// erdos_search3.cpp -- fast search for the six equation types (2.15)-(2.20)
//   v0 + v1 = v2 with each v either z^e (pure power, base not a perfect power)
//   or b_e(z)=z(z+1)...(z+e-1), 1/k+1/m+1/n < 1.   Same output as erdos_search2.
//
// Algorithm: for every exponent triple, the two SPARSEST of the three value sets are
// enumerated and the third (densest) set is queried by integer root extraction.
// Cost ~ N^(1/e1+1/e2) <= N^(7/12) instead of N.
//
// Build:  g++ -O3 -march=native -fopenmp -o erdos_search3 erdos_search3.cpp
// Usage:  ./erdos_search3 N [type=all|2.15|...|2.20] [--nogcd]      (summands <= N <= 4e18)

#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
typedef unsigned long long u64;
typedef unsigned __int128 u128;

static const u128 SAT = (u128)1 << 100;
static const int MODS[] = {64, 63, 65, 67, 71, 73};
static const int NM = 6;
static const int EB = 20;
static bool RAW = false;                         // reduced scan: perfect-power bases allowed
static int RK=0,RP=0,RQ=0;                        // fixed triple in reduced scan                       // max block length (21! > 4e18)

struct Fam { bool blk; int e; };                // pure with e==0 : the value 1 ("1^M")

static inline u128 evalF(const Fam& f, u64 z) {
    u128 r = 1;
    for (int i = 0; i < f.e; i++) {
        u64 t = f.blk ? z + i : z;
        if (r > SAT / t) return SAT;
        r *= t;
    }
    return r;
}

static bool isPP(u64 x) {                       // is x >= 2 a perfect power?
    for (int e = 2; e < 64; e++) {
        if (((u128)1 << e) > x) break;
        long double r = powl((long double)x, 1.0L / e);
        long long c = (long long)(r + 0.5L);
        for (long long d = max(2LL, c - 1); d <= c + 1; d++) {
            u128 v = 1; bool over = false;
            for (int i = 0; i < e; i++) { v *= (u64)d; if (v > x) { over = true; break; } }
            if (!over && v == x) return true;
        }
    }
    return false;
}

// number of values <= limit (approximate for pure powers: perfect-power bases ignored)
static u64 countVals(Fam f, u64 limit) {
    if (!f.blk && f.e == 0) return 1;
    u64 lo = 0, hi = 4000000000ULL;             // largest x with f(x) <= limit
    while (lo < hi) {
        u64 mid = (lo + hi + 1) / 2;
        if (evalF(f, mid) <= limit) lo = mid; else hi = mid - 1;
    }
    return lo;
}

struct Gen {
    vector<u64> val, base;
    vector<uint8_t> res[NM];
};
static map<tuple<bool,int,u64>, Gen*> cacheG;
static Gen* getGen(Fam f, u64 limit) {
    auto key = make_tuple(f.blk, f.e, limit);
    auto it = cacheG.find(key);
    if (it != cacheG.end()) return it->second;
    Gen* g = new Gen();
    if (!f.blk && f.e == 0) { g->val.push_back(1); g->base.push_back(1); }
    else for (u64 x = f.blk ? 1 : 2;; x++) {
        u128 v = evalF(f, x);
        if (v > limit) break;
        if (!RAW && !f.blk && isPP(x)) continue;
        g->val.push_back((u64)v); g->base.push_back(x);
    }
    for (int mi = 0; mi < NM; mi++) {
        g->res[mi].resize(g->val.size());
        for (size_t i = 0; i < g->val.size(); i++) g->res[mi][i] = (uint8_t)(g->val[i] % MODS[mi]);
    }
    return cacheG[key] = g;
}

struct Mem {                                    // membership test for one family
    Fam f; u64 fact = 1;
    vector<vector<char>> ok; int ord[NM];
    void init(Fam ff) {
        f = ff; fact = 1;
        if (f.blk) for (int i = 2; i <= f.e; i++) fact *= i;
        ok.assign(NM, {}); double dens[NM];
        for (int mi = 0; mi < NM; mi++) {
            int M = MODS[mi]; ok[mi].assign(M, 0); int cnt = 0;
            if (!f.blk && f.e == 0) { ok[mi][1 % M] = 1; cnt = 1; }
            else for (int z = 0; z < M; z++) {
                u128 r = 1;
                for (int i = 0; i < f.e; i++) r = (r * (u64)(f.blk ? z + i : z)) % M;
                if (!ok[mi][(int)r]) { ok[mi][(int)r] = 1; cnt++; }
            }
            dens[mi] = (double)cnt / M; ord[mi] = mi;
        }
        sort(ord, ord + NM, [&](int a, int b) { return dens[a] < dens[b]; });
    }
    bool find(u64 S, u64& z) const {
        if (!f.blk && f.e == 0) { z = 1; return S == 1; }
        if (f.blk && S % fact) return false;
        long double r = powl((long double)S, 1.0L / f.e);
        long long lo = (long long)r - (f.blk ? f.e + 1 : 2);
        if (lo < 1) lo = 1;
        long long hi = (long long)r + 2;
        for (long long c = lo; c <= hi; c++) {
            u128 v = evalF(f, (u64)c);
            if (v == S) {
                if (!f.blk && (c < 2 || (!RAW && isPP((u64)c)))) return false;
                z = (u64)c; return true;
            }
            if (v > S) break;
        }
        return false;
    }
};

static bool hyp(int a, int b, int c) {
    long long num = 0, den = 1;
    for (int e : {a, b, c}) if (e) { num = num * e + den; den *= e; }
    return num < den;
}
static u64 factorial(int e) { u64 r = 1; for (int i = 2; i <= e; i++) r *= i; return r; }
struct Type { const char* name; bool a, b, c; };

static void printTerm(bool blk, int e, u64 base, ostream& os) {
    if (blk) os << "b_" << e << "(" << base << ")";
    else if (e == 0) os << "1^M";
    else os << base << "^" << e;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s N [type] [--nogcd]\n", argv[0]); return 1; }
    u64 N = (u64)strtold(argv[1], nullptr);
    if (N > 4000000000000000000ULL) { fprintf(stderr, "N too large\n"); return 1; }
    string only = "all"; bool nogcd = false;
    if (getenv("RED")) { RAW = true; sscanf(getenv("RED"), "%d,%d,%d", &RK, &RP, &RQ); only = "2.17"; }
    for (int i = 2; i < argc; i++) { string s = argv[i]; if (s == "--nogcd") nogcd = true; else only = s; }
    int LA = 1; while (((u128)1 << (LA + 1)) <= N) LA++;
    int LC = 1; while (((u128)1 << (LC + 1)) <= (u128)2 * N) LC++;
    const u64 LIM[3] = {N, N, 2 * N};

    const Type types[] = {
        {"2.15", false, false, false}, {"2.16", false, false, true},
        {"2.17", true,  false, false}, {"2.18", true,  false, true},
        {"2.19", true,  true,  false}, {"2.20", true,  true,  true}};

    for (const Type& T : types) {
        if (only != "all" && only != T.name) continue;
        const bool blk[3] = {T.a, T.b, T.c};
        auto exps = [&](int role) {
            vector<int> v;
            if (blk[role]) for (int e = 2; e <= EB; e++) v.push_back(e);
            else { if (role != 2) v.push_back(0); for (int e = 2; e <= (role == 2 ? LC : LA); e++) v.push_back(e); }
            return v;
        };
        vector<int> E0 = exps(0), E1 = exps(1), E2 = exps(2);
        bool symAB = (T.a == T.b);
        string tn = T.name;
        for (int k0 : E0) for (int k1 : E1) {
            if (symAB && k0 > k1) continue;
            for (int k2 : E2) {
                if (RAW) { if (k0 != RK || k1 != RP || k2 != RQ) continue; }
                else if (!hyp(k0, k1, k2)) continue;
                Fam fam[3] = {{T.a, k0}, {T.b, k1}, {T.c, k2}};
                if (T.c && factorial(k2) > (u64)2 * N) continue;
                u64 cnt[3];
                bool empty = false;
                for (int r = 0; r < 3; r++) {
                    cnt[r] = countVals(fam[r], LIM[r]);
                    if (cnt[r] == 0) empty = true;
                }
                if (empty) continue;
                int t = 0;                                   // target = densest set
                for (int r = 1; r < 3; r++) if (cnt[r] > cnt[t] || (cnt[r] == cnt[t] && r == 2)) t = r;
                int p = (t == 0) ? 1 : 0, q = (t == 2) ? 1 : 2;
                if (t == 1) { p = 0; q = 2; }
                Gen* GP = getGen(fam[p], LIM[p]);
                Gen* GQ = getGen(fam[q], LIM[q]);
                Mem mt; mt.init(fam[t]);
                u64 g0 = 1;
                if (!nogcd) {
                    if (tn == "2.18") g0 = min(factorial(k0), factorial(k2));
                    if (tn == "2.19" || tn == "2.20") g0 = min(factorial(k0), factorial(k1));
                }
                bool same = symAB && k0 == k1;
                long long nP = GP->val.size(); size_t nQ = GQ->val.size();
                #pragma omp parallel for schedule(dynamic, 4)
                for (long long i = 0; i < nP; i++) {
                    u64 vp = GP->val[i];
                    for (size_t j = 0; j < nQ; j++) {
                        u64 vq = GQ->val[j];
                        u64 vt;
                        // roles: v0 + v1 = v2
                        if (t == 2) vt = vp + vq;                    // p=0,q=1
                        else {                                       // p is v0 (t==1) or v1 (t==0); q = v2
                            if (vq <= vp) continue;
                            vt = vq - vp;
                            if (vt > N) continue;
                        }
                        bool good = true;
                        for (int s = 0; s < NM; s++) {
                            int mi = mt.ord[s], M = MODS[mi];
                            int a = GP->res[mi][i], b = GQ->res[mi][j], x;
                            if (t == 2) { x = a + b; if (x >= M) x -= M; }
                            else { x = b - a; if (x < 0) x += M; }
                            if (!mt.ok[mi][x]) { good = false; break; }
                        }
                        if (!good) continue;
                        u64 bt;
                        if (!mt.find(vt, bt)) continue;
                        u64 base[3], val[3];
                        base[p] = GP->base[i]; val[p] = vp;
                        base[q] = GQ->base[j]; val[q] = vq;
                        base[t] = bt;          val[t] = vt;
                        u64 A = val[0], B = val[1];
                        if (same && A > B) continue;
                        u64 g = std::gcd(A, B);
                        if (!nogcd) {
                            if (g != g0) continue;
                            if (tn == "2.17" && k0 == 4 && B == 1) continue;
                        }
                        ostringstream os;
                        os << T.name << " | ";
                        printTerm(T.a, k0, base[0], os); os << " + ";
                        printTerm(T.b, k1, base[1], os); os << " = ";
                        printTerm(T.c, k2, base[2], os);
                        os << " | " << A << " + " << B << " = " << val[2] << " | gcd=" << g << "\n";
                        #pragma omp critical
                        { cout << os.str(); cout.flush(); }
                    }
                }
            }
        }
    }
    return 0;
}
