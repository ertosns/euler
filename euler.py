'''
Author ertosns
Date 23/8/2017

euler challenge solutions.
'''
import math

def palindrome(int_p):
    str_p = str(int_p)
    l = len(str_p)
    for i in range(l/2):
        if str_p[i] != str_p[l-i-1]:
            return False
    return True

#todo enhance
def prime(n):
    if n < 1:
        return False
    elif n<=3 :
        return True
    elif not n%2 or not n%3:
        return False
    i = 5;
    while i <= math.sqrt(n): 
        if not n%(i-2) and not n%i and not n%(i+2):
            return False
        i+=6
    return True

def get_nth_prime(n):
    i = 1
    pcount = 0
    while True:
        if prime(i):
            pcount+=1
        if pcount == n:
            return i
        i+=1

#sieve of erathostenes
#Complexity: O(Nlog(N)log(N))
#memorty: O(N)
#optimization:use segmanted algorithms for better memory performance for large enough number.
def eratosthenes(N):
    primes = []
    cnt = 0
    sieve = [True for num in range(N+1)]
    trav = 2
    ind = 0
    while math.pow(trav,2) <= N: #traversal [2,3,4, ..
        inc = 0
        while True: #incremental 
            ind = int(math.pow(trav,2)+inc*trav)
            if (ind > N):
                break
            sieve[ind-1] = False
            inc+=1
        trav+=1    
    
    for num in sieve:
        cnt+=1
        if num:
            primes.append(cnt)
            
    return primes

def prime_sieve(N, algo):
    return algo(N)
        
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
    t = int(raw_input())
    for x0 in xrange(t):
        N = int(raw_input())
        palind = -1
        tota = N
        p_range = 999

        for p0 in range(p_range, 99, -1):
            for p1 in range(p_range, 99, -1):
                tota = p0*p1
                if tota>N:
                    continue
                elif (palindrome(tota)) and tota>palind:
                    palind = tota
            p_range-=1
        print palind
#-----------------#
euler4()
