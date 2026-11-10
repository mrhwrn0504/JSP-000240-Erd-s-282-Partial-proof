// INDEPENDENT certificate in (a, b mod M) coordinates, using the plain odd-greedy step on integers:
//   q = smallest odd integer >= b/a ; r = a q - b ; g = gcd(r, b q) ; a' = r/g ; b' = b q / g.
// Edge (a, b mod M) -> (a', b' mod M) iff some integer lift b realises it with 2<=a'<=A.
// b' mod M depends only on b mod M*a*r (and r on b mod 2a), so lifts mod lcm(M,2a)*a*(2a) suffice
// (we simply enumerate all residues of b modulo W = M*2a*a*R where R covers every possible r via lcm? -> instead two-stage).
#include <bits/stdc++.h>
using namespace std; typedef long long ll; typedef __int128 LL;
static LL gcdL(LL x,LL y){ if(x<0)x=-x; if(y<0)y=-y; while(y){LL t=x%y;x=y;y=t;} return x; }
int main(int argc,char**argv){ int A=atoi(argv[1]); ll M=atoll(argv[2]);
  // states: a in [2,A], c in [0,M)
  vector<vector<char>> alive(A+1,vector<char>(M,1)); ll total=0;
  for(int a=2;a<=A;a++)for(ll c=0;c<M;c++){ if(__gcd((ll)a,__gcd(c,M))!=1 && false){} total++; }
  int round=0;
  while(true){ ll killed=0;
    #pragma omp parallel for schedule(dynamic) reduction(+:killed)
    for(int a=2;a<=A;a++){
      ll M1=M/__gcd(M,(ll)2*a)*(2*a); // lcm(M,2a): b mod M1 fixes r
      for(ll c=0;c<M;c++){ if(!alive[a][c])continue; bool ok=false;
        for(ll c1=c;c1<M1&&!ok;c1+=M){ if(c1%2==0) continue; if(__gcd((ll)a,c1%a)!=1) continue; // b must be odd, coprime to a
          // determine r from a representative
          LL b0=(LL)c1+(LL)M1*4; LL q=(b0+a-1)/a; if(q%2==0)q++; ll r=(ll)(a*q-b0); if(r<=0||r>=2*a){fprintf(stderr,"bad r\n");exit(1);}
          LL W=(LL)M1*a*r; // b mod W fixes b' mod M
          for(LL b=b0%W; b<W && !ok; b+=M1){ LL bb=b+W; // positive representative, large
            LL qq=(bb+a-1)/a; if(qq%2==0)qq++; LL rr=a*qq-bb; if(rr!=r){fprintf(stderr,"r mismatch\n");exit(1);}
            LL g=gcdL(rr,bb*qq); ll an=(ll)(rr/g); if(an<2||an>A) continue; LL bn=bb/gcdL(bb,g); // careful: b q / g
            LL num=bb*qq/g; ll cn=(ll)(num%M); if(alive[an][cn]) ok=true; }
        }
        if(!ok){alive[a][c]=0;killed++;} } }
    total-=killed; round++; fprintf(stderr,"round %d killed %lld alive %lld\n",round,killed,total); if(killed==0||total==0)break; }
  printf("b-coordinates: A=%d M=%lld rounds=%d alive=%lld\n",A,M,round,total); }
