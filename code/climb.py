# Track the numerator sequence of the odd greedy orbit of a/b far beyond what exact integers allow,
# by carrying q_k modulo a shrinking smooth modulus (exact first state, then exact congruence dynamics):
#   g=gcd(r,q^2); a'=r/g; r'=(-a q^2/g) mod 2a'; q'=(a q^2 - r q + g r')/r   (mod M/gcd(M,r))
import sys
from math import gcd, factorial
def primerange(lo,hi):
    s=bytearray([1])*(hi+1)
    for i in range(2,hi):
        if s[i]:
            if i>=lo: yield i
            for j in range(i*i,hi+1,i): s[j]=0
def run(a,b,P,E,T):
    W=1
    for p in primerange(2,P+1):
        e=max(E, int(E*600/p)+2)
        W*=p**e
    q=(b+a-1)//a
    if q%2==0:q+=1
    r=a*q-b; M=W; q%=M
    nums=[a]
    for k in range(T):
        if r==0: return nums,'TERMINATED'
        if M%(2*r): return nums,'precision lost (r=%d)'%r
        g=gcd(r,q*q%r if q*q%r else r); an=r//g
        if an==1: nums.append(1); return nums,'TERMINATED'
        rn=(-a*((q*q//g)%(2*an)))%(2*an)
        # q*q//g exact needs q^2 divisible by g: q known mod M (r|M) so fine
        N=a*q*q-r*q+g*rn
        d=gcd(r,M); Mn=M//d
        q=(N//d)*pow(r//d,-1,Mn)%Mn
        M=Mn; a,r=an,rn; nums.append(a)
    return nums,'step limit'
if __name__=='__main__':
    P=int(sys.argv[1]);E=int(sys.argv[2]);T=int(sys.argv[3])
    for N in range(int(sys.argv[4]),int(sys.argv[5])+1):
        nums,st=run(3,factorial(N)-1,P,E,T)
        runs=[];cur=1
        for i in range(1,len(nums)):
            if nums[i]==nums[i-1]+1:cur+=1
            else:runs.append(cur);cur=1
        runs.append(cur)
        print(N,'steps',len(nums)-1,'max',max(nums),'longest +1 run',max(runs),st,'| tail',nums[-8:],flush=True)
