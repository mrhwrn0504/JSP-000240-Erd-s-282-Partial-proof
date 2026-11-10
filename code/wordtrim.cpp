// Trim for a fixed periodic numerator word (all g=1) modulo M.
// state (k mod n, q mod M); edge if some lift q+M*i (0<=i<r) satisfies gcd(q,r)=1, (-a q^2) mod 2r == r', and maps to alive state.
#include <bits/stdc++.h>
using namespace std; typedef long long ll; typedef __int128 LL;
int main(int argc,char**argv){
  ll M=atoll(argv[1]); vector<int> w; for(int i=2;i<argc;i++) w.push_back(atoi(argv[i]));
  int n=w.size(); ll H=M/2; vector<vector<char>> alive(n,vector<char>(H,1)); ll total=(ll)n*H; int round=0;
  while(true){ ll killed=0; auto nxt=alive;
    for(int k=0;k<n;k++){ int a=w[k],r=w[(k+1)%n],rn=w[(k+2)%n];
      for(ll j=0;j<H;j++){ if(!alive[k][j])continue; ll q0=2*j+1; bool ok=false;
        for(int i=0;i<r&&!ok;i++){ LL q=(LL)q0+(LL)M*i; LL q2=q*q; if(__gcd((ll)r,(ll)(q2%r))!=1)continue;
          ll t=(ll)(q2%(2*r)); ll rr=((-(ll)a*t)%(2*r)+2*r)%(2*r); if(rr!=rn)continue;
          LL N=(LL)a*q2-(LL)r*q+rn; LL qn=N/r; ll qm=(ll)(qn%M); if(alive[(k+1)%n][(qm-1)/2]) ok=true; }
        if(!ok){nxt[k][j]=0;killed++;} } }
    alive.swap(nxt); total-=killed; round++; if(killed==0||total==0)break; }
  printf("M=%lld rounds=%d alive=%lld of %lld\n",M,round,total,(ll)n*H);
}
