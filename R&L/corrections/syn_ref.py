import numpy as np
import matplotlib.pyplot as plt
import math

mu = 0.012277471

for i in [1, 2, 3, 4]:
	fig=plt.figure()

	ax=fig.add_subplot(2,2,1)

	t,x,y = np.loadtxt("correction"+str(i)+".dat", usecols=(0,1,2), unpack=True)
	ax.plot(x,y, linewidth='1')
	ax.scatter(-mu, 0, label='m1', color='gold')
	ax.scatter((1-mu), 0, label='m2', color='red')
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	ax.set_title("Synodic RF, method  "+r'$\delta x_0=0$'+", vector "+str(i))

	ax=fig.add_subplot(2,2,2)

	t,x,y = np.loadtxt("correction"+str(i)+".dat", usecols=(0,1,2), unpack=True)
	r = np.sqrt(x**2 + y**2)
	theta = np.arctan2(y,x)

	theta1 = theta + t

	x1 = r * np.cos(theta1)
	y1 = r * np.sin(theta1)

	ax.plot(x1,y1, linewidth='1')
	ax.plot(-mu*np.cos(t), -mu*np.sin(t), label='m1', color='gold')
	ax.plot((1-mu)*np.cos(t), (1-mu)*np.sin(t), label='m2', color='red')
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_title("Inertial RF, method "+r'$\delta x_0=0$'+", vector "+str(i))


	ax=fig.add_subplot(2,2,3)

	t,x,y = np.loadtxt("correction2"+str(i)+".dat", usecols=(0,1,2), unpack=True)
	ax.plot(x,y, linewidth='1')
	ax.scatter(-mu, 0, label='m1', color='gold')
	ax.scatter((1-mu), 0, label='m2', color='red')
	ax.set_xlabel("x")
	ax.set_ylabel("y")

	ax.set_title("method "+r'$\delta v_{y_0}=0$')

	ax=fig.add_subplot(2,2,4)
	t,x,y = np.loadtxt("correction2"+str(i)+".dat", usecols=(0,1,2), unpack=True)

	r = np.sqrt(x**2 + y**2)
	theta = np.arctan2(y,x)

	theta1 = theta + t

	x1 = r * np.cos(theta1)
	y1 = r * np.sin(theta1)

	ax.plot(x1,y1, linewidth='1')
	ax.plot(-mu*np.cos(t), -mu*np.sin(t), label='m1', color='gold')
	ax.plot((1-mu)*np.cos(t), (1-mu)*np.sin(t), label='m2', color='red')
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	#ax.legend(loc='center right', bbox_to_anchor=(1, 0.5))
	ax.set_title("method "+r'$\delta v_{y_0}=0$')


	plt.tight_layout()
	plt.savefig("syn_orb_"+str(i)+".png")
	plt.show()