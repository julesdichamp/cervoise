import itertools
import random
import fastrand

tab1 = [i for i in range(0, 11)]
tab2 = [j for j in range(11, 21)]

caca = list(itertools.combinations(tab2, 2))
print(caca)

caca1 = list(itertools.combinations(tab1, 1))
print(caca1)

#pipi = [list(zip(x,tab1)) for x in itertools.permutations(caca,len(tab1))]
pipi = [ x + y for x in caca1 for y in caca]
print(pipi)

a = {}
a[(2,5,7,8)] = 0
a[(1,0)] = 1
print("AAAAAAAAAAAA", len(a), sum(list(a.values())))

b = ()
print(b + (2,4) + (3,))


