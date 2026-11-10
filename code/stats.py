import random, sys
from collections import Counter
from math import gcd
from climb import run
random.seed(7); NS=int(sys.argv[1]); bits=int(sys.argv[2])
L=Counter(); mx=Counter(); lost=0; longest=(0,None)
for _ in range(NS):
    b=random.getrandbits(bits)|1|(1<<(bits-1)); a=random.randrange(2,200)
    if gcd(a,b)!=1: continue
    nums,st=run(a,b,3000,2,400)
    if st!='TERMINATED': lost+=1; continue
    n=len(nums)-1; L[n]+=1; mx[max(nums).bit_length()]+=1
    if n>longest[0]: longest=(n,(a,b,max(nums)))
tot=sum(L.values()); print('terminated',tot,'not resolved',lost,'longest',longest[0],'max numerator there',longest[1][2])
c=0
for n in sorted(L):
    c+=L[n]
    if n%3==0 or tot-c<20: print(n, 'P(length>%d)=%.5f'%(n,(tot-c)/tot))
