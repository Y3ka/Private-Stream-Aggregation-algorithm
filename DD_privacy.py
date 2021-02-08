# Print a graph to compare the distributed differential privacy scheme against a naive scheme where geometric noise is added to all data participants independantly
from PSA_charm_scheme import *

x = [10**i for i in range(1,3)] #values for x-axis
print(x)
y = [] #average error for x participants for our DD privacy scheme
y2= [] #average error for x participants for the naive scheme
result = []
nb_samples = 100 # average calculated over nb_samples runs
for k in x:
    nb_users = k
    mean = 0
    mean2 = 0
    for j in range(0, nb_samples):
        sk = keygen(nb_users,p)
        result = error(nb_users, sk, data_group)
        mean += result[0]/nb_samples
        mean2 += result[1]/nb_samples
    y2.append(mean2)
    y.append(mean)

print(f"x = {x}")
print(f"average error for our DD privacy scheme = {y}")
print(f"average error for naive scheme = {y2}")

fig, ax = plt.subplots()
ax.plot(x, y, linestyle='-', marker='o', color='b', label="Our scheme")
ax.plot(x, y2, linestyle='--', marker='o', color='r', label="Naive scheme")
ax.set_xlabel("Number of users")
ax.set_ylabel("Error")
ax.legend()
plt.show()
