# Print a graph to compare the distributed differential privacy scheme against a naive scheme where geometric noise is added to all data participants independantly
from PSA_charm_scheme import *
import matplotlib.pyplot as plt

x = [10**i for i in range(1,4)] #abscisse
y = [] #average error for x participants for our DD privacy scheme
y2= [] #average error for x participants for the naive scheme
y3 = [] #upper bound value
result = []
nb_samples = 50 # mean calculated over nb_samples runs
bound = 4*delta/eps*(sqrt((1/gamma)*log(1/d)*log(2/0.9)))
for k in x:
    nb_users = k
    mean = 0
    mean2 = 0
    for j in range(0, nb_samples):
        sk = keygen(nb_users,p)
        result = gen_error(nb_users, sk, data_group)
        mean += result[0]/nb_samples
        mean2 += result[1]/nb_samples
    y2.append(mean2)
    y.append(mean)
    y3.append(bound)

print(f"x = {x}")
print(f"average error for our DD privacy scheme = {y}")
print(f"average error for naive scheme = {y2}")
print(f"upper bound value for the error = {bound}")

fig, ax = plt.subplots()
ax.plot(x, y, linestyle='-', marker='o', color='b', label="Our scheme")
ax.plot(x, y2, linestyle='--', marker='o', color='r', label="Naive scheme")
ax.plot(x, y3, linestyle='-', marker='o', color='g', label="upper bound")
ax.set_xlabel("Number of users")
ax.set_ylabel("Error")
ax.legend()
plt.semilogx()
plt.show()
