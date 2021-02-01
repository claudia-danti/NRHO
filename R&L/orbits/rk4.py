import numpy as np
import matplotlib.pyplot as plt

mu=0.012277471
alpha=pow(((1./3.)*mu/(1.-mu)),1./3.)
gamma1 = alpha*(1.-(1./3.)*alpha-(1./9.)*pow(alpha,2)-(23./81.)*pow(alpha,3))
gamma2 = alpha*(1.+(1./3.)*alpha-(1./9.)*pow(alpha,2)-(31./81.)*pow(alpha,3))
gamma3 = 1.-(7./12.)*mu*(1.+(23./84.)*pow((7./12.)*mu,2))

print(1/2-mu, np.sqrt(3)/2, 'L4')
print(1/2-mu, -np.sqrt(3)/2, 'L5')
print(1-mu-gamma1,0, 'L1')
print(1-mu+gamma2,0, 'L2')
print(-mu-gamma3,0, 'L3')


fig=plt.figure()

ax = fig.add_subplot()
data = np.loadtxt("rk4.dat")
ax.plot(data[:,1],data[:,2])
ax.scatter(1/2-mu, np.sqrt(3)/2, label='L4', color='orange')
ax.scatter(1/2-mu, -np.sqrt(3)/2, label='L5',color='gold')
ax.scatter(1-mu-gamma1,0, label = 'L1')
ax.scatter(1-mu+gamma2,0, label = 'L2')
ax.scatter(-mu-gamma3,0, label = 'L3')
ax.scatter(-mu,0, label='Earth', color='black')
ax.scatter(1-mu,0, label='Moon', color='brown')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.legend()
ax.set_title("rk4 after one period orbit")

plt.tight_layout()
plt.savefig("rk4.png")
plt.show()
