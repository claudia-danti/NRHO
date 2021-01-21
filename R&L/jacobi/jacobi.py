import numpy as np
import matplotlib.pyplot as plt


fig=plt.figure()


ax=fig.add_subplot(2,2,1)
data = np.loadtxt("jacobi11.dat")
ax.plot(data[:,0],data[:,1], label='method 1')
data= np.loadtxt("jacobi21.dat")
ax.plot(data[:,0],data[:,1], label='method 2')

ax.set_ylabel("C")
ax.legend()
ax.set_title("Jacobi constant, vector 1")

ax=fig.add_subplot(2,2,2)
data = np.loadtxt("jacobi12.dat")
ax.plot(data[:,0],data[:,1], label='method 1')
data= np.loadtxt("jacobi22.dat")
ax.plot(data[:,0],data[:,1], label='method 2')
ax.legend()
ax.set_ylabel("C")
ax.set_title("Jacobi constant, vector 2")

ax=fig.add_subplot(2,2,3)
data = np.loadtxt("jacobi13.dat")
ax.plot(data[:,0],data[:,1], label='method 1')
data= np.loadtxt("jacobi23.dat")
ax.plot(data[:,0],data[:,1], label='method 2')
ax.set_xlabel("t")
ax.set_ylabel("C")
ax.legend()
ax.set_title("Jacobi constant, vector 3")

ax=fig.add_subplot(2,2,4)
data = np.loadtxt("jacobi14.dat")
ax.plot(data[:,0],data[:,1], label='method 1')
data= np.loadtxt("jacobi24.dat")
ax.plot(data[:,0],data[:,1], label='method 2')
ax.set_xlabel("t")
ax.set_ylabel("C")
ax.legend()
ax.set_title("Jacobi constant, vector 3")

#plt.tight_layout()
plt.savefig("jacobi1.png")
plt.show()
