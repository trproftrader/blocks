// gdiv.cpp: all solutions of  b_k(x) = U^g - V^g  (g>=3, V>=1, gcd(U,V)=1) for x<=X.
// Every solution has d=U-V | b_k(x); for each divisor d (d^g<b) V is found by bisection on the
// increasing function f(V)=sum_{j=1}^g C(g,j) d^j V^{g-j} = b   (g=3: square test).
// Build: g++ -O3 -march=native -fopenmp -o gdiv gdiv.cpp    Usage: ./gdiv k g X
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
typedef unsigned long long u64; typedef unsigned __int128 u128;
static const u128 CAP=((u128)1)<<126;
static u128 smul(u128 a,u128 b){ if(a==0||b==0) return 0; if(a>CAP/b) return CAP; return a*b; }
static u128 spow(u128 a,int e){ u128 r=1; for(int i=0;i<e;i++) r=smul(r,a); return r; }
static u128 binom[30][30];
static int G;
// f(V) = (V+d)^g - V^g computed as positive sum, early exit when > lim
static u128 fval(u128 V,u128 d,u128 lim){
  u128 s=0;
  for(int j=1;j<=G;j++){ u128 t=smul(smul(binom[G][j],spow(d,j)),spow(V,G-j)); s+=t; if(s>lim||s>=CAP) return CAP; }
  return s;
}
static u128 isqrt128(u128 n){ u128 s=(u128)sqrtl((long double)n); while(s*s>n) s--; while((s+1)*(s+1)<=n) s++; return s; }
static u64 iroot(u64 v,int e){ long double r=powl((long double)v,1.0L/e); u64 c=(u64)(r+0.5L);
  for(u64 t=(c>2?c-2:1);t<=c+2;t++){ u128 p=1; bool ov=false; for(int i=0;i<e;i++){p*=t; if(p>v){ov=true;break;}} if(!ov&&p==v) return t; } return 0; }
// v = base^e with e maximal
static void maxpow(u64 v,u64&base,int&e){ base=v;e=1; if(v<2) return; for(int ee=63;ee>=2;ee--){ if(((u128)1<<ee)>v) continue; u64 r=iroot(v,ee); if(r){base=r;e=ee;return;} } }
int main(int argc,char**argv){
  int K=atoi(argv[1]); G=atoi(argv[2]); u64 X=strtoull(argv[3],0,10);
  for(int n=0;n<30;n++){binom[n][0]=1;for(int j=1;j<=n;j++)binom[n][j]=binom[n-1][j-1]+(j<=n-1?binom[n-1][j]:0);}
  vector<uint32_t> spf(X+K+2); iota(spf.begin(),spf.end(),0);
  for(u64 i=2;i*i<=X+K+1;i++) if(spf[i]==i) for(u64 j=i*i;j<=X+K+1;j+=i) if(spf[j]==j) spf[j]=i;
  long long cnt=0;
  #pragma omp parallel for schedule(dynamic,2048)
  for(long long x=1;x<=(long long)X;x++){
    map<u64,int> f; u128 b=1;
    for(int t=0;t<K;t++){ u64 n=x+t; b*=n; while(n>1){u64 p=spf[n]; f[p]++; n/=p;} }
    vector<u128> ds={1};
    for(auto&pe:f){ size_t m=ds.size(); u128 pw=1; for(int i=1;i<=pe.second;i++){ pw*=pe.first; for(size_t j=0;j<m;j++) ds.push_back(ds[j]*pw);} }
    for(u128 d:ds){
      if(spow(d,G)>=b) continue;
      u128 V=0;
      if(G==3){ u128 q=b/d; u128 disc=12*q-3*d*d; if(12*q<=3*d*d) continue; u128 s=isqrt128(disc); if(s*s!=disc) continue; if((s-3*d)%6) continue; V=(s-3*d)/6; }
      else { long double hi_ld=powl((long double)(b/(G*d)),1.0L/(G-1)); u128 lo=1,hi=(u128)hi_ld+3; bool found=false;
        while(lo<=hi){ u128 mid=(lo+hi)/2; u128 fv=fval(mid,d,b); if(fv==b){V=mid;found=true;break;} if(fv<b) lo=mid+1; else hi=mid-1; }
        if(!found) continue; }
      if(V<1) continue;
      u128 U=V+d;
      if(fval(V,d,b)!=b) continue;
      u64 gg=std::gcd((u64)V,(u64)U); if(gg!=1) continue;
      u64 bv,bu; int ev,eu; maxpow((u64)V,bv,ev); maxpow((u64)U,bu,eu);
      #pragma omp critical
      { printf("k=%d g=%d x=%lld V=%llu U=%llu | V=%llu^%d U=%llu^%d | b=%.6e\n",K,G,x,(u64)V,(u64)U,(unsigned long long)bv,ev,(unsigned long long)bu,eu,(double)b); fflush(stdout); }
    }
  }
}
