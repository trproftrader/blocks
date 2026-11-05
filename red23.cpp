// red23.cpp : all solutions of  b_k(x) + v^2 = w^3  (x,v,w >= 1, w^3 <= Nmax)
// b_k(x) = x(x+1)...(x+k-1).  For each w, every x with b_k(x) < w^3 is tried; w^3 - b_k(x)
// must be a perfect square.
// Build: g++ -O3 -march=native -fopenmp -o red23 red23.cpp      Usage: ./red23 k Nmax
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
typedef unsigned long long u64;
typedef unsigned __int128 u128;

static bool sqmod[4][80];
static const int M[4] = {64, 63, 65, 61};

int main(int argc, char** argv) {
    int k = atoi(argv[1]);
    u64 N = (u64)strtold(argv[2], nullptr);
    for (int i = 0; i < 4; i++) for (int r = 0; r < M[i]; r++) sqmod[i][(u64)r * r % M[i]] = true;
    // all x with b_k(x) <= N
    vector<u64> bv;
    for (u64 x = 1;; x++) {
        u128 p = 1; for (int i = 0; i < k; i++) p *= (x + i);
        if (p > N) break;
        bv.push_back((u64)p);
    }
    u64 W = (u64)cbrtl((long double)N) + 2;
    while ((u128)W * W * W > N) W--;
    #pragma omp parallel for schedule(dynamic, 64)
    for (long long wi = 2; wi <= (long long)W; wi++) {
        u64 w = wi, w3 = w * w * w;
        for (size_t xi = 0; xi < bv.size(); xi++) {
            u64 b = bv[xi];
            if (b >= w3) break;
            u64 D = w3 - b;
            if (!sqmod[0][D % 64] || !sqmod[1][D % 63] || !sqmod[2][D % 65] || !sqmod[3][D % 61]) continue;
            u64 s = (u64)sqrtl((long double)D);
            while ((u128)s * s > D) s--;
            while ((u128)(s + 1) * (s + 1) <= D) s++;
            if ((u128)s * s == D) {
                #pragma omp critical
                { printf("%d %zu %llu %llu %llu\n", k, xi + 1, (unsigned long long)s, (unsigned long long)w, (unsigned long long)b); fflush(stdout); }
            }
        }
    }
    return 0;
}
