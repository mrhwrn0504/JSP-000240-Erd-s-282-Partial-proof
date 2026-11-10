#!/usr/bin/env python3
"""Per-prime FORALL-trim for a periodic numerator word (all g=1).
State (k, q mod p^s).  A state is kept iff EVERY lift q+p^s*i satisfies the step-k condition
   a q^2 + r' = 0 (mod p^{v_p(2r)})   [a=w_k, r=w_{k+1}, r'=w_{k+2}]
and its image q'=(a q^2 - r q + r')/r (mod p^s) is kept at phase k+1.
A non-empty greatest fixed point = a residue class ALL of whose integers follow the word forever at p."""
import sys
from padic import vp, primes_of
def forall_trim(word,p,s):
    n=len(word); P=p**s
    qs=[q for q in range(P) if (p!=2 or q%2==1)]
    alive=[set(qs) for _ in range(n)]
    # precompute images
    img=[{} for _ in range(n)]
    for k in range(n):
        a=word[k]; r=word[(k+1)%n]; rn=word[(k+2)%n]
        e=vp(2*r,p); er=vp(r,p)
        for q in qs:
            res=set(); ok=True
            for i in range(p**er):
                Q=q+P*i
                if (a*Q*Q+rn)%p**e: ok=False;break
                N=a*Q*Q-r*Q+rn
                qn=(N//p**er)*pow(r//p**er,-1,P)%P
                res.add(qn)
            img[k][q]=res if ok else None
    changed=True
    while changed:
        changed=False
        for k in range(n):
            nx=alive[(k+1)%n]
            dead=[q for q in alive[k] if img[k][q] is None or not img[k][q]<=nx]
            if dead:
                changed=True
                alive[k].difference_update(dead)
    return alive
if __name__=='__main__':
    word=tuple(map(int,sys.argv[1].split(',')))
    for p in primes_of(word):
        for s in range(1,int(sys.argv[2])+1):
            al=forall_trim(word,p,s)
            print(p,s,[len(x) for x in al])
