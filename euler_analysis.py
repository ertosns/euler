import euler
import math
import time
'''
algorithmic dynamic analysis for euler algorithms
'''
def kthfreepower_complexity_wrapper(li):
    periods = []
    for k in range(2, int(math.log(2, li))):
        #the use of global segmanted sieve would corrrupt the calculations.
        prime_init = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
        start = time.time()
        res = euler.kthfreepower(li, k)
        stop = time.time()
        periods+=[stop-start]
        #verify res.
    return euler.avg(periods)
def kthfreepower_complexity():
    L = euler.random_list(1000000, 10)
    P = []
    O = []
    for l in L:
        if l<2:
            continue
        P += [kthfreepower_complexity_wrapper(l)]
    for i in range(len(L)):
        O+=[math.log(float(P[i+1])/P[i], float(L[i+1])/L[i])]
        i+=1
    return euler.avg(O)


#-----------------#
#TODO impl grenade dogwatch.
power_limit = 4
#time.time() fails and return difference of 0
assert(kthfreepower_complexity()<power_limit)
