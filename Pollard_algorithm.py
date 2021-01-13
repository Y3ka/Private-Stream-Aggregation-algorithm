#Implementation of Pollard algorithm for our charm crypto scheme
from charm.toolbox.integergroup import integer
from charm.core.math.integer import reduce

def naturals_from(n):
    while True:
        yield n
        n += 1           
                
def step_xab(x, a, b, alpha, beta, n, Z):
    s = integer(x, 3)

    # S1
    if s == integer(1, 3):
        x = reduce(integer(integer(x) * integer(beta), Z))
        b = reduce(integer(integer(b) + integer(1), n))
        
    
    # S2
    if s == integer(0, 3):
        x = reduce(integer(integer(x) * integer(x), Z))
        a = reduce(integer(integer(2) * integer(a), n))
        b = reduce(integer(integer(2) * integer(b), n))

    # S3
    if s == integer(2, 3):
        x = reduce(integer(integer(x) * integer(alpha), Z))
        a = reduce(integer(integer(a) + integer(1), n))
    
    return x, a, b
    

def pollard_rho(alpha, beta, n, Z):
    x = {0: integer(1)}
    a = {0: integer(0)}
    b = {0: integer(0)}
    
    # returns x, a, b for a given i using memoization
    def get_xab(i):
        if i not in x:
            _x, _a, _b = get_xab(i - 1)
            x[i], a[i], b[i] = step_xab(_x, _a, _b, alpha, beta, n, Z)
        return x[i], a[i], b[i]
        
    #print("i\tx_i\ta_i\tb_i\tx_2i\ta_2i\tb_2i")
    for i in naturals_from(1):
        x_i, a_i, b_i = get_xab(i)
        x_2i, a_2i, b_2i = get_xab(2 * i)
        
        #print("%d\t%d\t%d\t%d\t%d\t%d\t%d" % (i, x_i, a_i, b_i, x_2i, a_2i, b_2i))
        if x_i == x_2i:
            r = reduce(integer(integer(b_i)-integer(b_2i), n))
            if r == 0:
                return False
            else:
                return reduce(integer(integer(1/r) * integer(reduce(integer(integer(a_2i) - integer(a_i), n))), n))