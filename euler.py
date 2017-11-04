'''
Author ertosns
Date 23/8/2017

euler challenge solutions.

guidlines:
- algorithms could coplicated for the sake of performance.

uncompleted:
- 7, fails the last two  cases, after trying eratosthenes sieve, and current solution with different variations,

'''
import math
import time

prime_init = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
def mult(nums):
    mul = 1
    for i in nums:
        mul*=i
    return mul

def avg(l):
    if not len(l):
        return 0
    i = 0
    sum = 0
    while i < len(l):
        sum+=l[i++]
    return float(sum)/i

#time analysis, for several modes.
#hardlimite|watchdog, analyis relation between speed, and process resources, it's application in limited systems like android.
#TODO create regression test for computational complexity, given  complexity as equation, and input validate implementation.
def palindrome(int_p):
    str_p = str(int_p)
    l = len(str_p)
    for i in range(l/2):
        if str_p[i] != str_p[l-i-1]:
            return False
    return True

def prime_sieve(N, algo):
    pcount = 1
    prime = 1
    primes = []
    for index in algo([0, N]):
        if index:
            primes+=[prime]
            pcount+=1
        prime+=1             
    return primes

def prime(n):
    if n < 1:
        return False
    elif n<=3 :
        return True
    elif not n%2 or not n%3:
        return False
    i = 5;
    while i <= math.sqrt(n):
        # 5+6k (3+6k|3) -> <5+6k, 7+6k>
        if  not n%i or not n%(i+2):
            return False
        i+=6
    return True

# illuminate non-primes by nonlinear combinations of initial prime list.
# time complexity O(N)
# below non working example of two multiplicatives, how achieve arbitrary level?
def ertosns(N):
    P = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    p0 = P[0]
    plen = len(P)
    plast = P[len(P)-1]
    plimit = p0*plast+1
    S = [True]*N
    for i in range(1, int(N/2)+1):
        S[i*2-1] = False
    bound = 0
    pb = int(math.log(plast, p0))
    O_pb = 1
    index = plast+1
    res = 1
    pw = 1
    compl = True
    
    while compl:
        shift = bound
        for Pi in P[:bound+1]: #old
            for power in range(O_pb, pb+1):
                pw = int(math.pow(Pi, power))
                for Pj in P[shift+1:]:
                    res = pw*Pj
                    if res<N:
                        S[res-1] = False
                if pw<N and power>1: #primes powers
                    S[pw-1] = False
            shift+=1
        shift = bound
        for Pi in P[bound:]: #new
            for power in range(1, pb+1):
                pw = int(math.pow(Pi, power))
                for Pj in P[shift+1:]:
                    if power==1 and Pi*Pj>N:
                        compl = False
                    res = pw*Pj
                    if res<N:
                        S[res-1] = False
                if pw<N and power>1: #new powers
                    S[pw-1] = False
            shift+=1
            
        
        for pot in S[plast: plimit]:
            if pot:
                P+=[index]
            index+=1
        
        bound = plen-1
        oplen = plen
        plen = len(P)
        plast = P[plen-1]
        index = plast+1
        plimit = p0*plast+1
        if plimit>=N:
            plimit = N-1
        O_pb=pb
        pb = (int(math.log(plast, p0)))
        
    return [1, 2]+P

def ertosns_test():
    primes = ertosns(100000)
    for p in primes:
        if not prime(p):
            print 'not a prime '+str(p)
            break

#sieve of erathostenes
#Complexity: O(Nlog(N)log(N))
#memorty: O(N)
def eratosthenes(N):
    pass

def eratosthenes_test():
    pass


''' 
complexity: ?
l (list of primes) should be passed as parameter for more than one use.
'''

#todo update new name usage
def pth_trial(nth, l):
    if nth<len(l):
        return l
    sus = l[len(l)-1]
    while len(l)<nth:
        sus+=2 #insert wheel sieving.
        floor = math.floor(sus**0.5)
        for p in l:
            if sus%p == 0:
                break
            if p>floor:
                l.append(sus)
                break
    return l

def ntrial(N, l):
    last = l[len(l)-1]
    if N<last:
        return l
    while last<N:
        last+=2 #insert wheel sieving.
        floor = math.floor(last**0.5)
        for p in l: #limit p_range
            if last%p == 0:
                break
            if p>floor:
                l.append(last)
                break
    return l

def kthfreepower(N, k):
    cfree = N
    for n in range (2, N+1):
        lim = n**(float(1)/k)
        ntrial(lim, prime_sieve)
        for p in prime_init:
            if not n%int(math.pow(p, k)):
                cfree-=1
                break
    return cfree

#not the use of global segmanted sieve in kthfreepower would corrrupt the calculations.
def kthfreepower_test(li):
    periods = []
    for k in range(2, int(math.log(2, li))):
        prime_init = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
        start = time.time()
        res = kthfreepower(li, k)
        stop = time.time()
        periods+=[stop-start]
        #verify res.
    return avg(periods)

def sieve_test():
    primes = trial(100000, prime_init)
    for p in primes:
        if not prime(p):
            print('not a prime '+str(p))
            break

'''
#TODO replace primes with stream of bits
#still fails the last two cases!
def get_nth_prime(n):
    primes = [2,3]
    cprime = len(primes)
    primes+=[0 for i in range(n)]
    prime = 5
    ptwo=True
    lprime = len(primes)
    while cprime<n:
        sqrt = math.sqrt(prime)
        for p in primes:
            if not prime%p:
                break
            if p>=sqrt:
                primes[cprime]=prime
                cprime+=1
                break
        
        prime+= 2 if ptwo else 4 #illuminate 4/6 of numbers, multiples of 2,3.
        ptwo = False if ptwo else True
        if cprime==lprime:
            primes+=[0 for i in range(math.sqrt(n))]
            lprime = len(primes)
    return primes[n-1]
'''

#fermat algorithm
def fermat_fact(N):
    fact = []
    if N < 1:
        return fact
    if not N%2:
        return [2, N/2]
    q = int(math.ceil(math.sqrt(N)))
    while q<N:
        p2 = q*q-N
        p = int(math.sqrt(p2))
        if p2/p==p and (q+p) != N:
            fact.append([(q+p),(q-p)])
        q+=1
    return fact

#assumes N>1, return prime factors for N
def trial_division(N):
    fact = []
    for p in prime_sieve(N, eratosthenes):
        if not N%p:
            fact.append(p)
    return fact

def largest_prime_factor(N):
    n=N
    d = 2
    while n>1 and d<=math.sqrt(N):
        if not (n%d):
            n=n/d
        else:
            d = d+2 if d!=2 else 3
    if n==1:
        return d
    else:
        return n
    
def factorize(N, algo):
    return algo(N)

def has_ndig_factors(N, w):
    UL = math.pow(10,w)
    LL = math.pow(10,w-1)-1
    for pair in factorize(N, fermat_fact):
        if pair[0] < UL and pair[0] > LL and pair[1] < UL and pair[1] > LL:
            return True
    return False

#--testing--#
def test():
    #eratosthenes_test()
    #ertosns_test()def avg(l):
    i = 0
    sum = 0
    while i < len(l):
        sum+=l[i++]
    return float(sum)/i
    sieve_test()

#todo grenade dogwatch.
def complexity():
    L = [1, 2, 3] #updated with varied lengths
    P = []
    O = []
    for l in L:
        P += [kthfreepower_test(l)]
    for i in range(len(L)):
        O+=[math.log(float(P[i+1])/P[i], float(L[i+1])/L[i])]
        i+=1
    return avg(O)
    
#--analysis--#

#TOOD how receieve, doted-args in python?
def computational_time(algorithm, inpts):
    #fill profile, and run regression, for function coefficients.
    pass

#--problem panel--#

#find the sum of all numbers below N multiple of 3,5.
def euler1():
    t = int(raw_input())
    for a0 in range(t):
        n = int(raw_input())
        i = 1
        sum = 0
        while i*3 < n:
            sum += i*3
            if not i%3 and i*5<n:
                sum += i*5
            i+=1
        print (sum)
#sum even fabonacci under N
def euler2():
    t = int(raw_input())
    for a0 in xrange(t):
        n=long(raw_input())
        if n<2:
            print 0
        elif n==2:
            print 2
        else:
            i=1
            j=2
            sum=0
            while i<n:
                if i%2==0:
                    sum+=i
                tmp=j
                j=i+j
                i=tmp
            print sum
#max prime factor under N
def euler3():
    t = int(raw_input())
    for x0 in xrange(t):
        n = int(raw_input())
        print largest_prime_factor(n)
#search max palindrome number under N as product of 3-digits
def euler4():
    '''
    1  2  3 ...  N
    | /| /|     /|
    |/ |/ |    / |
    1  2  3 ...  N
    '''
    t = int(raw_input())
    for x0 in xrange(t):
        N = int(raw_input())
        found = 0
        for i in range(999, 99, -1):
            for j in range(100, i+1):
                palind = i*j
                if palind>=N or palind < found:
                    continue
                if palindrome(palind):
                    found = palind
                    
        print found
#smallest multiple of all numbers under N
def euler5():
    t = int(raw_input())
    for a0 in range(t):
        N = int(raw_input())
        primes = prime_sieve(N, eratosthenes)
        mult = 1
        for p in primes:
            if p>1:
                x = int(math.log(N,p))
                mult *= int(math.pow(p,x))
        print mult
#difference between square of sums, and sums of squares
def euler6():
    t = int(raw_input())
    for a0 in range(t):
        N = int(raw_input())
        print int(math.pow(N*(N+1)/2, 2))-(N*(N+1)*(2*N+1)/6)

#get nth prime
#reusability of generated primes is substantial.
def euler7():
    t = int(input())
    for a0 in range(t):
        N = int(input())
        if N>len(l):
            trial(N, prime_init)
        print prime_init[N-1]

def euler8():
    t = int(input())
    for a0 in range(t):
        line0 = raw_input().split()
        N = raw_input()
        limits = [int(i) for i in line0]
        digs = []
        for i in range(limits[0]):
            digs+=[int(N[i])]
            
        poss = []
        for i in range(len(digs)-limits[1]+1):
            mult = 1
            tmp = digs[i:limits[1]+i]
            for dig in tmp:
                mult*=dig
            poss+=[mult]
        poss.sort()
        print (poss[len(poss)-1])
'''
given N=a+b+c "sum of pythagorean triplet" validate that,
and return  a*b*c, -1 otherwise. 
(1) N=a+b+c
(2) c*c=a*a+b*b
----------------
(3) b=(N*N-2*N*a)/(2*N-2*a)
'''
def euler9():
    t = int(input())
    for a0 in range(t):
        N = int(input())
        a = 1
        product=-1
        while a<=N/3:
            b = int((N*N - 2*N*a)/(2*N-2*a))
            c = N-a-b
            if c*c==b*b+a*a:
                product=a*b*c
                break
            a+=1
        print(product)
#just like euler7, but reusability of old work Sums solves the problem.
def euler10():
    t = int(input())
    prime_init = [3, 5, 7, 11, 13]
    Sum = [0, 0, 2, 5, 5, 10, 10, 17, 17, 17, 17, 28, 28] #under 13
    slen = 0
    OP = 1
    N = 1
    for a0 in range(t):
        N = int(input())
        OP = len(prime_init)-1
        trial(int(float(N)/math.log(N)), prime_init)
        while prime_init[len(prime_init)-1]<=N:
            trial(len(prime_init)*10, prime_init)
        if len(Sum) <= prime_init[len(prime_init)-1]:
            while (OP+1)<len(prime_init):
                Sum += [Sum[len(Sum)-1]+prime_init[OP]]* \
                       (prime_init[OP+1]-prime_init[OP])
                OP+=1
        print Sum[N]

def euler11():
    mtrx = []
    for a0x in range(20):
        mtrx += [[int(i) for i in raw_input().split()]]
    trav = 0
    mult = 0
    rc = 0
    r = 0
    while r<20:
        c = 0
        while c<20:
            rc = mtrx[r][c]
            if c<=16: #forward
                trav = rc*mtrx[r][c+1]*mtrx[r][c+2]*mtrx[r][c+3]
                if trav>mult:
                    mult = trav
            if r<=16: #down
                trav = rc*mtrx[r+1][c]*mtrx[r+2][c]*mtrx[r+3][c]
                if trav>mult:
                    mult = trav
            if r<=16 and c<=16: #right diagonal
                trav = rc*mtrx[r+1][c+1]*mtrx[r+2][c+2]*mtrx[r+3][c+3]
                if trav>mult:
                    mult = trav
            if r<=16 and c>=3: #left diagonal
                trav = rc*mtrx[r+1][c-1]*mtrx[r+2][c-2]*mtrx[r+3][c-3]
                if trav>mult:
                    mult=trav
            c+=1
        r+=1
    print mult

def euler12():
    t0 = int(input())
    for a0 in range(t0):
        n = int(input())
        primes = primes(n)
        start = mult(primes)
        N = int((math.sqrt(1-8*start)-1)/2) #positive quadratic formula
        tr=N
        while True:
            div = divisors(tr) #todo, implement, return all divisors.
            if len(div)>n:
                print tr
                break
            tr+=N
            N+=1

def euler13():
    t0 = int(input())
    sum =0
    for a0 in range(t0):
        sum+=int(input())
    print str(sum)[0:10]
    
#number of primes is huge, so segmented sieve has to be used.
def euler193():
    args = [int(i) for i in raw_input().split()]
    N = args[0]
    pf = N
    K = args[1]    
    print kthfreepower(N, K)

#-----------------#
#test()
#complexity()
#euler193()
ntrial(int(input()), prime_init)

#update problems with complexity.
#write recommendations
#write list of used unproven theorems
#use partition set count (bell number) fact in driving all combination of factors.
#factors
#radical factors
#prime factors
#ertosns above is uncomplete sieving by illumination of nonprimes by nonlinear combination of initial prime list, how calculate nonlinear combination?
