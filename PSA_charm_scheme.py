#Implementation for the Private Stream Aggregation algorithm proposed by Elaine Shi
#Based on www.elaineshi.com/docs/ndss2011.pdf
from charm.toolbox.integergroup import IntegerGroup
from charm.toolbox.integergroup import integer
from charm.toolbox.hash_module import Hash
from charm.toolbox.conversion import Conversion
from math import exp, log,sqrt
from datetime import datetime
import hashlib, secrets
import scipy.stats as st
from Pollard_algorithm import *


def H(t):
    """
        Compute hash value for time step t (in G)
    """
    global G
    return G.hash(t)



#We subclass the generic discrete random variable class
class Geometric(st.rv_discrete):
    """
        Two-sided geometric distribution
        Probability mass function at k, l > 1: f(k) = (l-1)/(l+1).l**(-|k|)
        One sided : f(k) = (l-1).l**k
    """
    def _pmf(self, k, l):
        # return (l-1)*pow(l, -k)
        return (((l-1)/(l+1))*pow(l, -abs(k)))



class Aggregator():
    """
        Computing and decrypting the sum of users' data
    """
    def __init__(self, ciphertexts, sk):
        #list of encrypted value for time t
        self.c = ciphertexts
        #aggregator's secret key
        self.sk = sk

    def decrypt(self):
        """
            Decrypt the sum.
            Return the decrypted sum for time step t
        """
        global g, ht, delta
        v1 = pow(ht,self.sk)
        v2 = prod_enc(self.c)
        v = (v1*v2)
        #sum_dec = compute_log(g, v, delta)
        sum_dec = integer(pollard_rho(g, v, p, modulus))
        return sum_dec

def prod_enc(c):
    """
        Compute the product of all users' ciphertexts for time step t
    """
    prod = c[0]
    for i in range(1,len(c)):
        prod *= c[i]
    return prod

def compute_log(g, v, delta):
    """
        Brute-force search for computing the discrete log
    """
    global nb_users
    for sum in range (0, (nb_users*delta)+1):
        value = pow(g, sum)
        if v == value:
            return sum
    print("Error during discrete log computation: no value found")
    exit()

class User():
    """
        Computing the noise and encrypting user's data
    """
    def __init__(self, userdata, sk, nb_users):
        global g, ht, eps, gamma, d, delta
        #users' data
        self.x = userdata 
        #user' secret key
        self.sk = sk
        #hash value
        self.ht = ht 
        #generator
        self.g = g 
        #number of users
        self.nb_users = nb_users
        self.alpha = exp(eps/delta)
        self.beta = (1/(gamma*nb_users))*log(1/d)

    def encrypt(self):
        """
            Compute user's ciphertext for time step t (NoisyEnc(param,ski,t,x))
            input :
            g : public param - prime order, generator
            sk : user's secret key
            ht : hash value for time t
            x : value to encrypt
            return :
            c : user's ciphertext for time t
        """
        c1 = pow(ht,self.sk)
        c2 = pow(g,self.x)
        c = (c1*c2)
        return c

    def DD_Privacy(self):
        """
            Implement the DD-Private Data Randomization Procedure
            Generating and adding the noise necessary for achieving distributed differential privacy
        """
        bernoulli_dist = st.bernoulli(self.beta) #Bernouilli's law with parameter beta
        dist = Geometric(a=-self.x, b=delta-self.x+1, name="Geometric")
        if bernoulli_dist.rvs():
            r = dist.rvs(l=self.alpha)
        else:
            r = 0
        self.x = (self.x + r)%p
    
    def naive_scheme(self):
        """
            Each user add independant geometric noise to their data
        """
        dist = Geometric(a=-self.x, b=delta-self.x+1, name="Geometric")
    
        r = dist.rvs(l=self.alpha)
        self.x = (self.x + r)%p

def keygen(n, p):
    """
        Choose n+1 random secrets [sk0, sk1,..., skn] such that sk0+sk1+...+skn = 0 mod p
        return :
        L : list of the secret keys
    """
    L = []
    sum = 0
    for i in range(0, n):
        number  = 0
        if sum < -p:
            l = [i for i in range(-p - sum, p+1)]
            number = secrets.choice(l)
        elif sum > p:
            l = [i for i in range(-p, p - sum+1)]
            number = secrets.choice(l)
        else:
            l = [i for i in range(-p , p+1)]
            number = secrets.choice(l)
        sum += number
        L.append((number+p)%p)
    number = (-sum)%p
    L.append(number)

    return L

def gen_error(nb_users, sk, data_group):
    """
        Compute the error created by the DD_Privacy and the error if using a naive geometric noise
        return a list with these errors: [error, error_naive]
    """
    error = 0
    error_naive = 0
    users = []
    users_naive = []
    for i in range(1, nb_users+1):
        user = User(secrets.choice(data_group), sk[i], nb_users)
        users.append(user)
        users_naive.append(user)

    #list with the clear values
    data_list = []

    #Users add noise to their data
    for i in range(0, len(users)):
        data_list.append(users[i].x)
        users[i].DD_Privacy()
        error += (users[i].x - data_list[i])
        users_naive[i].naive_scheme()
        error_naive += (users_naive[i].x - data_list[i])
        
    return [error, error_naive]

#############################################
# TESTS AND SIMULATIONS
#############################################

#group definition (G = Schnorr group)
G = IntegerGroup()
#number of bits used for prime numbers p and q
bits = 14
#Generate p and q
G.paramgen(bits)

#prime order
p = int(G.q)
print(f"\nGroup order: {p}")
modulus = 2*p+1
print(f"Modulus: {modulus}")

#generator
g = G.randomGen()
print(f"Group generator: {g}")

#Hash value for time t
t = 1
ht = H(t)

#text space
#delta : limit value of user's data (user's value x[i] in {0,1,...,delta})
delta = 20
data_group = [i for i in range(0,delta+1)]

# (eps, d)-DD-privacy parameters with delta >= eps/3 and 0 < d < 1
eps = 0.5
d = 0.1

#proportion of uncompromised participants
gamma = 1


#number of users
nb_users = 10
print(f"Number of users: {nb_users}")

#secret keys
sk = keygen(nb_users,p)
sum = 0
for i in sk:
    sum = (sum+i)%p


print(f"secret keys : {sk}")
print(f"sum of the secret keys : {sum}\n")

#############################################
# Test for a given time t and n users
#############################################

#Creation of users
users = []
for i in range(1, nb_users+1):
    user = User(secrets.choice(data_group), sk[i], nb_users)
    users.append(user)

data_list = []
expected_result = 0
expected_result_with_noise = 0
ciphertexts = []
ciphertexts_noisy = []

#Users add noise to their data
for user in users:
    data_list.append(user.x)
    expected_result = (expected_result + user.x)%p
    ciphertexts.append(user.encrypt())
    user.DD_Privacy()
    expected_result_with_noise = (expected_result_with_noise + user.x)%p
    ciphertexts_noisy.append(user.encrypt())

#Aggregator
agg = Aggregator(ciphertexts, sk[0])
agg_noisy = Aggregator(ciphertexts_noisy, sk[0])
decrypt = agg.decrypt()
decrypt_noisy = agg_noisy.decrypt()
error = abs(int(decrypt) - int(decrypt_noisy))
print(f"list of all users'data (before encryption) : {data_list}")
print(f"list of all ciphertexts (given to the aggregator) : {ciphertexts}\n")
print(f"expected sum : {expected_result}")
# print(f"expected sum (with noise): {expected_result_with_noise}")
print(f"decrypted sum by the aggregator : {decrypt}")
# print(f"decrypted sum by the aggregator with noise : {decrypt_noisy}")
# print(f"Error range : {error}")
