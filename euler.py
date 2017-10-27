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

def mult(nums):
    mul = 1
    for i in nums:
        mul*=i
    return mul

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

#sieve of erathostenes
#Complexity: O(Nlog(N)log(N))
#memorty: O(N)
def eratosthenes(N):
    if N[1]<1:
        return []
    rng = N[1]-N[0]+1
    size = (rng)//2
    sieve = [True]*size
    limit = int(N[1]**0.5)
    for i in range(N[0]+1, limit):
        if sieve[i]:
            val = 2*i+1
            span = (size-1-i)//val
            sieve[i+val-N[0]::val] = [False]*span
    return sieve

def eratosthenes_test():
    for i in range(1000):
        i+=2
        lb=int(math.sqrt(i))         
        sieve = eratosthenes([lb,i])
        for j in range (i-lb):
            if prime(j+lb+1):
                assert(sieve[j])

def get_nth_prime(n):
    primes = []
    nsqrt = int(n**0.5)
    ON = 0
    N = n+nsqrt
    sieve = eratosthenes([ON, N])
    pcount = 1
    prime = 1
    while pcount<=n:
        if sieve[prime-ON]:
            primes+=[prime]
            pcount+=1
        prime+=2
        if prime-ON >= len(sieve)-1:
            ON = N
            N = N+nsqrt
            sieve = eratosthenes([ON, N])
    return primes[n-1]


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
    eratosthenes_test()

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
def euler7():
    t = int(input())
    for a0 in range(t):
        N = int(input())
        print get_nth_prime(N)

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
                    
def euler10():
    t = int(input())
    for a0 in range(t):
        N = int(input())
        primes = prime_sieve(N, eratosthenes)
        print sum(primes)

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

# complexity O(N**2/lnN)
def euler193():
    args = [int(i) for i in raw_input().split()]
    N = args[0]
    K = args[1]
    primes = prime_sieve(N, eratosthenes)[1:]
    sf = [True]*N #i*p**k == j*p`**k --> (i/j)**(1/k)*p = p`
    pk = 0
    for p in primes:
        pk = p**k
        i = 1
        while i*pk<=N:
            sf[i*pk] = False
            i+=1

    sfcount = 0
    for num in sf:
        if num:
            sfcount+=1
    print sfcount

#-----------------#
#test()
euler193()

#update problems with complexity.
#write recommendations
#write list of used unproven theorems
#factors
#radical factors
#prime factors
