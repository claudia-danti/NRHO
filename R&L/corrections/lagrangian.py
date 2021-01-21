import numpy as np
import matplotlib.pyplot as plt

mu=0.012277471
alpha=pow(((1./3.)*mu/(1.-mu)),1./3.)
gamma1 = alpha*(1.-(1./3.)*alpha-(1./9.)*pow(alpha,2)-(23./81.)*pow(alpha,3))
gamma2 = alpha*(1.+(1./3.)*alpha-(1./9.)*pow(alpha,2)-(31./81.)*pow(alpha,3))
gamma3 = 1.-(7./12.)*mu*(1.+(23./84.)*pow((7./12.)*mu,2))

c1 = [2.96795205692014, 2.70666498827617, 3.05718591486646, 3.07768490992014]
c2 = [2.98505025238624, 2.70805466542589, 3.05783274738714, 3.07992172496128]




for cont in [1, 2, 3, 4]:
	fig=plt.figure()

	ax=fig.add_subplot(1,2,1)
	x,y = np.loadtxt("correction"+str(cont)+".dat", usecols=(1,2), unpack=True)
	ax.plot(x,y)
	ax.scatter(1/2-mu, np.sqrt(3)/2, label='L4', color='orange')
	ax.scatter(1/2-mu, -np.sqrt(3)/2, label='L5',color='gold')
	ax.scatter(1-mu-gamma1,0, label = 'L1')
	ax.scatter(1-mu+gamma2,0, label = 'L2')
	ax.scatter(-mu-gamma3,0, label = 'L3')
	ax.scatter(-mu,0, label='m1', color='black')
	ax.scatter(1-mu,0, label='m2', color='brown')

	x_vals = np.linspace(-1.5, 1.5, 1000)
	y_vals = np.linspace(-1.5, 1.5, 1000)
	X, Y = np.meshgrid(x_vals, y_vals)

	Z=pow(X,2) + pow(Y,2) + 2*(1-mu)/np.sqrt(pow(X+mu,2)+pow(Y,2)) + 2*mu/np.sqrt(pow(X-1+mu,2)+pow(Y,2)) - c1[cont-1]

	ax.contour(X,Y,Z,0, colors = 'violet')
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.legend()
	ax.set_title("Corrected initial conditions, "+r'$\delta x_0=0$')

	ax=fig.add_subplot(1,2,2)
	data = np.loadtxt("correction2"+str(cont)+".dat")
	ax.plot(data[:,1],data[:,2])
	ax.scatter(1/2-mu, np.sqrt(3)/2, label='L4', color ='orange')
	ax.scatter(1/2-mu, -np.sqrt(3)/2, label='L5', color='gold')
	ax.scatter(1-mu-gamma1,0, label = 'L1')
	ax.scatter(1-mu+gamma2,0, label = 'L2')
	ax.scatter(-mu-gamma3,0, label = 'L3')
	ax.scatter(-mu,0, label='m1', color='black')
	ax.scatter(1-mu,0, label='m2', color='brown')
	x_vals = np.linspace(-1.5, 1.5, 1000)
	y_vals = np.linspace(-1.5, 1.5, 1000)
	X, Y = np.meshgrid(x_vals, y_vals)

	Z=pow(X,2) + pow(Y,2) + 2*(1-mu)/np.sqrt(pow(X+mu,2)+pow(Y,2)) + 2*mu/np.sqrt(pow(X-1+mu,2)+pow(Y,2)) - c2[cont-1]

	ax.contour(X,Y,Z,0, colors ='violet')
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.legend()
	ax.set_title("Corrected initial conditions, "+r'$\delta v_{y_0}=0$')


	plt.savefig("correction1.png")
	plt.show()

