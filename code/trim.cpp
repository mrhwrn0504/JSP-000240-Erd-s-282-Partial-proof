// Over-approximation graph of the odd-greedy dynamics mod M.
// State: (letter (a,r), q mod M, q odd).  Exact map on integers:
//   g=gcd(r,q^2); a'=r/g; r'=(-a q^2/g) mod 2a'; q'=(a q^2 - r q + g r')/r.
// Edge (a,r,q mod M)->(a',r',q' mod M) iff some lift q+M*i (0<=i<r) realises it with 2<=a'<=A.
// If iterated trimming kills every state, NO orbit with numerators in [2,A] of length > #rounds exists.
#include <bits/stdc++.h>
using namespace std;
typedef long long ll; typedef __int128 LL;
int A; ll M;
vector<pair<int,int>> letters; map<pair<int,int>,int> lid;
int main(int argc,char**argv){
  A=atoi(argv[1]); M=atoll(argv[2]);
  if(M%2){fprintf(stderr,"M must be even\n");return 1;}
  for(int a=2;a<=A;a++)for(int r=1;r<2*a;r++) if(__gcd(a,r)==1 && (a+r)%2==1){lid[{a,r}]=letters.size();letters.push_back({a,r});}
  int L=letters.size(); ll H=M/2; // odd residues q=2j+1
  vector<vector<bool>> alive(L, vector<bool>(H,true));
  ll total=(ll)L*H; int round=0; vector<pair<int,ll>> lastk,prevk;
  while(true){
    ll killed=0; vector<vector<bool>> nxt=alive; prevk.swap(lastk); lastk.clear();
    #pragma omp parallel for schedule(dynamic) reduction(+:killed)
    for(int l=0;l<L;l++){
      int a=letters[l].first, r=letters[l].second;
      for(ll j=0;j<H;j++){
        if(!alive[l][j]) continue;
        ll q0=2*j+1; bool ok=false;
        for(int i=0;i<r && !ok;i++){
          LL q=(LL)q0+(LL)M*i; LL q2=q*q;
          ll g=__gcd((ll)r,(ll)(q2% r)); // gcd(r,q^2)
          int an=r/g; if(an<2||an>A) continue;
          ll t=(ll)((q2/g)%(2*an)); ll rn=((-(ll)a*t)%(2*an)+2*an)%(2*an);
          if(rn==0) continue; // would terminate
          LL N=(LL)a*q2-(LL)r*q+(LL)g*rn; // divisible by r
          LL qn=N/r; ll qm=(ll)(qn%M);
          auto it=lid.find({an,(int)rn}); if(it==lid.end()) continue;
          if(alive[it->second][(qm-1)/2]) ok=true;
        }
        if(!ok){nxt[l][j]=false;killed++;
          #pragma omp critical
          { if((int)lastk.size()<200000) lastk.push_back({l,j}); } }
      }
    }
    alive.swap(nxt); total-=killed; round++;
    fprintf(stderr,"round %d killed %lld alive %lld\n",round,killed,total);
    if(killed==0||total==0)break;
  }
  printf("A=%d M=%lld rounds=%d alive=%lld\n",A,M,round,total);
  if(total==0 && getenv("DUMPLAST")){ FILE*f=fopen(getenv("DUMPLAST"),"w"); for(auto&x:lastk) fprintf(f,"%d %d %lld %d\n",letters[x.first].first,letters[x.first].second,2*x.second+1,round); for(auto&x:prevk) fprintf(f,"%d %d %lld %d\n",letters[x.first].first,letters[x.first].second,2*x.second+1,round-1); fclose(f);} 
  if(total>0 && getenv("ZSTAT")){ // histogram of (p, v_p(r), v_p(z)) over alive states, z=2aq-r
    int ps[]={2,3,5,7,11,13,17,19,23,29,31}; map<array<int,3>,ll> Hh;
    for(int l=0;l<L;l++){int a=letters[l].first,r=letters[l].second; for(ll j=0;j<H;j++){ if(!alive[l][j])continue; ll q0=2*j+1;
      for(int p:ps){ if(M%p) continue; int cap=0; ll m=M; while(m%p==0){m/=p;cap++;} if(p==2)cap++; ll z=2LL*a*q0-r; int vz=0; ll pw=1; for(int i=0;i<cap;i++)pw*=p; ll zz=((z%pw)+pw)%pw; if(zz==0) vz=99; else while(zz%p==0){zz/=p;vz++;} int vr=0; int rr=r; while(rr%p==0){rr/=p;vr++;} Hh[{p,vr,vz}]++; } } }
    for(auto&e:Hh) printf("  p=%d v(r)=%d v(z)=%d : %lld\n",e.first[0],e.first[1],e.first[2],e.second);
  }
  if(total>0 && argc>3){ // dump alive->alive letter edges with multiplicity
    map<array<int,4>,ll> E;
    for(int l=0;l<L;l++){int a=letters[l].first,r=letters[l].second;
      for(ll j=0;j<H;j++){ if(!alive[l][j])continue; ll q0=2*j+1;
        for(int i=0;i<r;i++){ LL q=(LL)q0+(LL)M*i; LL q2=q*q; ll g=__gcd((ll)r,(ll)(q2%r)); int an=r/g; if(an<2||an>A)continue;
          ll t=(ll)((q2/g)%(2*an)); ll rn=((-(ll)a*t)%(2*an)+2*an)%(2*an); if(rn==0)continue;
          LL N=(LL)a*q2-(LL)r*q+(LL)g*rn; LL qn=N/r; ll qm=(ll)(qn%M); auto it=lid.find({an,(int)rn}); if(it==lid.end())continue;
          if(alive[it->second][(qm-1)/2]) E[{a,r,an,(int)rn}]++; }}}
    FILE*f=fopen(argv[3],"w"); for(auto&e:E)fprintf(f,"%d %d %d %d %lld\n",e.first[0],e.first[1],e.first[2],e.first[3],e.second); fclose(f);
  }
  if(total>0){ // print letter histogram of survivors
    for(int l=0;l<L;l++){ll c=0;for(ll j=0;j<H;j++)c+=alive[l][j]; if(c)printf("  (%d,%d): %lld\n",letters[l].first,letters[l].second,c);}
  }
}
