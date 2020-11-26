#Implementation for of the Private Stream Aggregation algorithm proposed by Elaine Shi
#TODO implement Pollard's method

from math import exp, log
from datetime import datetime
import hashlib, math, secrets
import scipy.stats as st
# secrets is used for generating cryptographically strong random numbers
# import numpy.random as rd


#We subclass the generic discrete random variable class
class Geometric(st.rv_discrete):
    """
        Symmetric geometric distribution
        Probability mass function at k, l > 1: f(k) = (l-1)/(l+1).l**(-|k|)
        One sided : f(k) = (l-1).l**k
    """
    def _pmf(self, k, l):
        return (l-1)*pow(l, -k)
        # return (((l-1)/(l+1))*pow(l, -abs(k))) #(l**(-abs(k)))

def generators(n):
    """
        Return the first generator of Z/Zn
    """
    s = set(range(1, n))
    results = []
    for a in s:
        g = set()
        for x in s:
            g.add((a**x) % n)
        if g == s:
            return a
    return -1


def H(p):
    """
        Compute hash value for time t
    """
    t = datetime.now()
    #converts into string
    s = t.strftime("%B/%D/%Y/%H:%M:%S")
    #convert into bytes for the hash function
    e = s.encode()
    #hash data 
    ht = int.from_bytes(hashlib.sha256(e).digest(), byteorder='little')
    
    return ht%p

def NoisyEnc(p, g, sk, ht, x):
    """
        Compute user's ciphertext for time step t (NoisyEnc(param,ski,t,x))
        input :
        p, g : public param - prime order, generator
        sk : user's secret key
        ht : hash value for time t
        x : encrypted value
        return :
        c : user's ciphertext for time t
    """
    c1 = pow(ht,sk,p)
    c2 = pow(g,x,p)
    c = (c1*c2)%p
    return c

def decrypt(p, g, sk, ht, c):
    """
        Decrypt the sum
        inputs :
        p, g : public param - prime order, generator
        sk : aggregator's secret key
        ht : hash value for time t
        c : list of encrypted value for time t
        return :
        sum_dec : decrypted sum for time t
    """
    global delta
    v1 = pow(ht,sk,p)
    v2 = prod_enc(p, c)
    v = (v1*v2)%p
    sum_dec = compute_log(g, v, p, delta)
    return sum_dec

def prod_enc(p, c):
    """
        Compute the product of all users' ciphertexts for time t
    """
    prod = 1
    for i in c:
        prod = (prod*i)%p
    return prod

def compute_log(g, v, p, delta):
    global nb_user
    for sum in range (0, (nb_user*delta)+1):
        value = pow(g, sum)%p
        if v == value:
            return sum
    return -1

def keygen(n, p):
    """
        Choose n+1 random secrets [sk0, sk1,..., skn] such that sk0+sk1+...+skn = 0 mod(p-1)
        return :
        L : list of the secret keys
    """
    L = []
    q = p-1
    sum = 0
    for i in range(0, n):
        number  = 0
        if sum < -q:
            l = [i for i in range(-q - sum, q+1)]
            number = secrets.choice(l)
            #number = rd.randint(-q - sum, q)
        elif sum > q:
            l = [i for i in range(-q, q - sum)]
            number = secrets.choice(l)
            #number = rd.randint(-q, q - sum)
        else:
            l = [i for i in range(-q , q+1)]
            number = secrets.choice(l)
            #number = rd.randint(-q, q)
        sum += number
        L.append((number+q)%q)
    number = (-sum)%q
    L.append(number)

    return L


def DD_Privacy(eps, gamma, d, delta, x):
    """
        Implement the DD-Private Data Randomization Procedure
        Generation and adding of the noise necessary for achieving distributed differential privacy
        inputs : 
        eps, d : (eps, d)-DD-privacy parameters
        gamma : proportion of uncompromised participants (number of honest participants = gamma*len(x))
        delta : limit value of user's data (x[i] in {1,2,...,K})
        x : list of all participants' data in a certain period of time
        return :
        x_noise : list of users' data with noise in that period of time
    """
    n = len(x) #number of participants
    alpha = exp(eps/delta)
    beta = (1/(gamma*n))*log(1/d)
    noise = []
    x_noise = []
    bernoulli_dist = st.bernoulli(beta) #Bernouilli's law with parameter beta
    for i in range (0,n):
        if bernoulli_dist.rvs():
            dist = Geometric(a=1, b=1000, name="Geometric")
            r = dist.rvs(l=alpha)
        else:
            r = 0
        noise.append(r)
        x_noise.append((x[i]+r)%p)
    print("noise : ", noise)
    return x_noise


#############################################
# TEST
#############################################

#prime number
p = 10067

#number of user
nb_user = 20

#text space
delta = 20
data_group = [i for i in range(0,delta+1)]

# user data
x= []
for i in range(0, nb_user):
    x.append(secrets.choice(data_group))
# x = [10,7,10,7,10,12,1,4,3,10]
print(f"\noriginal data of {nb_user} users : {x}")

# (eps, d)-DD-privacy parameters
# delta >= eps/3
# 0 < d < 1
eps = 20
d = 0.01

#p roportion of uncompromised participants
gamma = 0.9

x_noisy = DD_Privacy(eps, gamma, d, delta, x)
print(f"noisy data of {nb_user} users : {x_noisy}\n")



#secret keys
sk = keygen(nb_user,p)
sum = 0
for i in sk:
    sum = (sum+i)%(p-1)

print("secret keys : ", sk)
print("sum of the secret keys : ", sum)

#choose a generator
g = generators(p)
if not g:
#    print(f"genrators of Z/Z{p}: {gens}")
    print("No generator")
    exit()


#Hash value for time t
ht = H(p)


x_enc = []
x_enc_noisy = []
for i in range(0,nb_user):
    x_enc.append(NoisyEnc(p, g, sk[i+1], ht, x[i]))
    x_enc_noisy.append(NoisyEnc(p, g, sk[i+1], ht, x_noisy[i]))
print(f"encrypted data received by the aggregator : {x_enc}")
print(f"encrypted data received by the aggregator (with noise) : {x_enc}\n")
result = decrypt(p, g, sk[0], ht, x_enc)
result_noisy = decrypt(p, g, sk[0], ht, x_enc_noisy)
print(f"decrypted sum by the aggregator : {result}")
print(f"decrypted sum by the aggregator (with noise) : {result_noisy}\n")

#Compute the true result to compare with the one from decrypted sum
true_result = 0
true_result_noisy = 0
for i in x:
    true_result += i%(p-1)
for i in x_noisy:
    true_result_noisy += i%(p-1)

print(f"Sum of the unencrypted data : {true_result}")
print(f"Sum of the unencrypted noisy data : {true_result_noisy}")


 