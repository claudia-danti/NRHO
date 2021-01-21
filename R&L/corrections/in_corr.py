import numpy as np
import matplotlib.pyplot as plt

mu=0.012277471
alpha=pow(((1./3.)*mu/(1.-mu)),1./3.)
gamma1 = alpha*(1.-(1./3.)*alpha-(1./9.)*pow(alpha,2)-(23./81.)*pow(alpha,3))
gamma2 = alpha*(1.+(1./3.)*alpha-(1./9.)*pow(alpha,2)-(31./81.)*pow(alpha,3))
gamma3 = 1.-(7./12.)*mu*(1.+(23./84.)*pow((7./12.)*mu,2))


fig=plt.figure()

ax=fig.add_subplot(1,3,1)
data = np.loadtxt("initial4.dat")
ax.plot(data[:,1],data[:,2])
ax.scatter(1/2-mu, np.sqrt(3)/2, label='L4', color='orange')
ax.scatter(1/2-mu, -np.sqrt(3)/2, label='L5',color='gold')
ax.scatter(1-mu-gamma1,0, label = 'L1')
ax.scatter(1-mu+gamma2,0, label = 'L2')
ax.scatter(-mu-gamma3,0, label = 'L3')
ax.scatter(-mu,0, label='m1', color='black')
ax.scatter(1-mu,0, label='m2', color='brown')

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Initial conditions")

ax=fig.add_subplot(1,3,2)
data = np.loadtxt("correction4.dat")
ax.plot(data[:,1],data[:,2])
ax.scatter(1/2-mu, np.sqrt(3)/2, label='L4', color='orange')
ax.scatter(1/2-mu, -np.sqrt(3)/2, label='L5',color='gold')
ax.scatter(1-mu-gamma1,0, label = 'L1')
ax.scatter(1-mu+gamma2,0, label = 'L2')
ax.scatter(-mu-gamma3,0, label = 'L3')
ax.scatter(-mu,0, label='m1', color='black')
ax.scatter(1-mu,0, label='m2', color='brown')
ax.set_xlabel("x")

ax.set_title("Corrected initial conditions, "+r'$\delta x_0=0$')

ax=fig.add_subplot(1,3,3)
data = np.loadtxt("correction24.dat")
ax.plot(data[:,1],data[:,2])
ax.scatter(1/2-mu, np.sqrt(3)/2, label='L4', color='orange')
ax.scatter(1/2-mu, -np.sqrt(3)/2, label='L5',color='gold')
ax.scatter(1-mu-gamma1,0, label = 'L1')
ax.scatter(1-mu+gamma2,0, label = 'L2')
ax.scatter(-mu-gamma3,0, label = 'L3')
ax.scatter(-mu,0, label='m1', color='black')
ax.scatter(1-mu,0, label='m2', color='brown')
ax.set_xlabel("x")

ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_title("Corrected initial conditions, "+r'$\delta v_{y_0}=0$')


plt.savefig("correction1.png")
plt.show()
