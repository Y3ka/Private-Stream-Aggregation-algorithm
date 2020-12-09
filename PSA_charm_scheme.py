from charm.toolbox.integergroup import IntegerGroup
from charm.toolbox.hash_module import Hash
from charm.toolbox.conversion import Conversion
from math import exp, log
from datetime import datetime
import hashlib, secrets
import scipy.stats as st

# G = IntegerGroup()
# G.paramgen(8)
# p = G.random()
# g = G.randomGen()
# t = 5
# order = G.groupOrder()

# print(f"{p}.{H(t)} = {p*H(t)}")

# print(order)



def H(t):
    """
        Compute hash value for time step t (in G)
    """
    global G
    return G.hash(t)



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
        sum_dec = compute_log(g, v, delta)
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
    def __init__(self, userdata, sk):
        global g, ht, nb_users, eps, gamma, d, delta
        #users' data
        self.x = userdata 
        #user' secret key
        self.sk = sk
        #hash value
        self.ht = ht 
        #generator
        self.g = g 
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
            return :
            x_noise : user's data with noise
        """
        bernoulli_dist = st.bernoulli(self.beta) #Bernouilli's law with parameter beta
        if bernoulli_dist.rvs():
            dist = Geometric(a=1, b=1000, name="Geometric")
            r = dist.rvs(l=self.alpha)
        else:
            r = 0
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
            # if l == []:
            #     number = 0
            number = secrets.choice(l)
        elif sum > p:
            l = [i for i in range(-p, p - sum)]
            # if l == []:
            #     number = 0
            number = secrets.choice(l)
        else:
            l = [i for i in range(-p , p+1)]
            # if l == []:
            #     number = 0
            number = secrets.choice(l)
        sum += number
        L.append((number+p)%p)
    number = (-sum)%p
    L.append(number)

    return L


#############################################
# TEST
#############################################

#group definition
G = IntegerGroup()
G.paramgen(64)

#prime order
p = int(G.q)
print(f"\nGroup order: {p}")


#generator
g = G.randomGen()
print(f"Group generator: {g}")

#number of users
nb_users = 30
print(f"Number of users: {nb_users}")

#Hash value for time t
t = 1
ht = H(t)

#text space
#delta : limit value of user's data (x[i] in {0,2,...,delta})
delta = 50
data_group = [i for i in range(0,delta+1)]

# (eps, d)-DD-privacy parameters with delta >= eps/3 and 0 < d < 1
eps = 20
d = 0.01

#proportion of uncompromised participants
gamma = 0.9

#secret keys
sk = keygen(nb_users,p)
sum = 0
for i in sk:
    sum = (sum+i)%p

print(f"secret keys : {sk}")
print(f"sum of the secret keys : {sum}\n")


#Creation of users
users = []
for i in range(1, nb_users+1):
    user = User(secrets.choice(data_group), sk[i])
    users.append(user)

data_list = []
expected_result = 0
ciphertexts = []
#Users add noise to their data
for user in users:
    # user.DD_Privacy()
    data_list.append(user.x)
    expected_result = (expected_result + user.x)%p
    ciphertexts.append(user.encrypt())

#Aggregator
agg = Aggregator(ciphertexts, sk[0])
print(f"list of all users'data (before encryption) : {data_list}")
print(f"list of all ciphertexts (given to the aggregator) : {ciphertexts}\n")
print(f"decrypted sum by the aggregator : {agg.decrypt()}")
print(f"expected sum : {expected_result}")