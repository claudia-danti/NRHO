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

for cont in [1, 2, 3, 4]:
	ax = fig.add_subplot(2,2,cont)
	data = np.loadtxt("correction"+str(cont)+".dat")
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
	ax.legend()
	ax.set_title("Corrected initial vector "+str(cont)+", dx0=0")

plt.tight_layout()
plt.savefig("corrected_dx0=0.png")
plt.show()

fig=plt.figure()

for cont in [1, 2, 3, 4]:
	ax = fig.add_subplot(2,2,cont)
	data = np.loadtxt("correction2"+str(cont)+".dat")
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
	ax.legend()
	ax.set_title("Corrected initial vector "+str(cont)+", dvy0=0")

plt.tight_layout()
plt.savefig("corrected_dvy0=0.png")
plt.show()

fig=plt.figure()

for cont in [1, 2, 3, 4]:
	ax = fig.add_subplot(2,2,cont)
	data = np.loadtxt("correction3"+str(cont)+".dat")
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

	ax.set_title("Corrected initial vector "+str(cont)+", dx(T/2)=0")

fig=plt.figure()

for cont in [1, 2, 3, 4]:
	ax = fig.add_subplot(2,2,cont)
	data = np.loadtxt("correction4"+str(cont)+".dat")
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

	ax.set_title("Corrected initial vector "+str(cont)+", dvy(T/2)=0")

plt.tight_layout()
plt.savefig("corrected_dvyT2)=0.png")
plt.show()