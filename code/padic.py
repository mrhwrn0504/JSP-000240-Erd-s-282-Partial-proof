#!/usr/bin/env python3
"""Per-prime realisability of a periodic numerator word with all g=1.
Condition at step k (a,r,r')=(w_k,w_{k+1},w_{k+2}):  a q^2 + r' = 0 mod p^{v_p(2r)}  (and p∤q if p|r),
update q' = (a q^2 - r q + r')/r  in Z_p.  DFS over p-adic digits of q_0."""
import sys
def vp(n,p):
    e=0
    while n%p==0: n//=p;e+=1
    return e
def sim_p(word,p,q,n,K):
    """q known mod p^n. returns ('dead',k) | ('need',k) | ('ok',)"""
    L=len(word)
    for k in range(K):
        a=word[k%L]; r=word[(k+1)%L]; rn=word[(k+2)%L]
        e=vp(2*r,p)
        if e==0:
            # no condition; q' = (a q^2 - r q + rn)/r mod p^n
            if n==0: continue
            q=(a*q*q-r*q+rn)*pow(r,-1,p**n)%p**n
            continue
        if n<e: return ('need',k)
        if (a*q*q+rn)%p**e: return ('dead',k)
        er=vp(r,p)
        N=a*q*q-r*q+rn       # known mod p^n ; divisible by p^er
        n2=n-er
        q=(N//p**er)*pow(r//p**er,-1,p**n2)%p**n2 if n2>0 else 0
        n=n2
    return ('ok',)
def realis_p(word,p,K,limit=300000):
    start=[(1,1)] if p==2 else [(c,1) for c in range(p)]
    stack=start;nodes=0;best=0;cnt=0
    while stack:
        q,n=stack.pop();nodes+=1
        if nodes>limit: return ('limit',best,nodes)
        res=sim_p(word,p,q,n,K)
        if res[0]=='ok': return ('ok',q,n,nodes)
        if res[0]=='dead': best=max(best,res[1]);continue
        for i in range(p): stack.append((q+i*p**n,n+1))
    return ('none',best,nodes)
def primes_of(word):
    ps={2}
    for a in word:
        m=a;d=2
        while d*d<=m:
            if m%d==0:
                ps.add(d)
                while m%d==0:m//=d
            d+=1
        if m>1:ps.add(m)
    return sorted(ps)
if __name__=='__main__':
    word=tuple(map(int,sys.argv[2].split(',')));K=int(sys.argv[1])
    for p in primes_of(word):
        print(p,realis_p(word,p,K))
