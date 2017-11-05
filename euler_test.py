import euler

def mult_test():
    print('mult testing')
    assert(euler.mult([2,3,4,109341238479745, 0, -1923874192])==0)
    assert(euler.mult([-10928374, -29837443987, -19873219874])<0)
    print('mult test done')
def avg_test():
    print('avg testing')
    assert(euler.avg([1,2,3,4])==float(10)/4)
    assert(euler.avg([0,-1,1,-2,2,4, 0, 0, 0])==float(4)/9)
    print('avg test done!')
def palindrome_test():
    print('palindrome testing')
    assert(euler.palindrome("12321"))
    assert(euler.palindrome("1221"))
    assert( not euler.palindrome("001001"))
    assert(euler.palindrome("0"))
    assert(euler.palindrome("1111"))
    print('palindrome test done!')
def prime_sieve_test():
    print('prime_sieve testing')
    for i in range(1000):
        sieve  = euler.prime_sieve(i, euler.eratosthenes)
        for p in sieve:
            assert(euler.prime(p))
    print('prime_sieve test done!')
def prime_test():
    print('prime testing')
    for n in range(1000):
        for d in range(2, int(n**0.5)+1):
            if not n%d:
                assert(not euler.prime(n))
                break
            if d==(n**0.5):
                assert(euler.prime(n))
    print('print test done!')
def ertosns_test():
    primes = euler.ertosns(100000)
    for p in primes:
        if not euler.prime(p):
            print 'not a prime '+str(p)
            break
def eratosthenes_test():
    print('eratosthenes testing')
    for i in range(1000):
        for p in  euler.prime_sieve(i, euler.eratosthenes):
            assert(euler.prime(p))
    print('eratosthenes test done!')
def pth_trial_test():
    print('pth_trial testing')
    primes = euler.pth_trial(1000, euler.prime_init)
    i = 0
    while i < len(primes)-1:
        assert(euler.prime(primes[i]))
        for np in range(primes[i]+1, primes[i+1]):
            assert(not euler.prime(np))
        i+=1
    print('pth_trial test done!')
    
#--testing--#
mult_test()
avg_test()
palindrome_test()
prime_sieve_test()
prime_test()
#ertosns_test() uncomplete yet 
eratosthenes_test()
pth_trial_test()
