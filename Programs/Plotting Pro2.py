from numpy import array, zeros, linspace, polyfit, arange, where, amin, reshape
import matplotlib.pyplot as plt
from Py_Functions import readarrays


#Error analysis of non-interacting electrons

X = array([4,5,7,10,15])
N = array([100,200,300,400,500])

lambdas = array([3,7,11,15])



Errors = [zeros(shape=(5,5)) for i in range(4)]
minindexes = zeros(shape=(4,2))

for k in range(4):
	minimum = 1000
	m = 0 ; n = 0	

	for i in range(5):		
		for j in range(5): 
			eig = readarrays("Quantum_Eigenpairs_N_%d_X_%d.txt" % (N[j],X[i]))[0][0][k]
			Errors[k][i][j] = abs(lambdas[k]-eig)

			if Errors[k][i][j] < minimum:
				m = i;n = j
				minimum = Errors[k][i][j]
	minindexes[k][0] = m ; minindexes[k][1] = n;
	print "%d , %d" % (m,n)
	print Errors[k][m][n]

print minindexes
print Errors

tablefile = open("errortable_1particle.txt","w+")

tablefile.write("Ideal values for N and Xmax:\n")
for i in range(4):
	 tablefile.write("Lambda %d: Rho_max = %d | N = %d\n" % (i+1,X[int(minindexes[i][0])],N[int(minindexes[i][1])])) 
tablefile.close()				

errortable = open("errortable complete.txt","w+")
errortable.write("Complete Errortable:\n\n")
for k in range(4):
	errortable.write("Errors for lambda %d = %d:\n" % (k+1,lambdas[k]))
	errortable.write("               N = 100     N = 200     N = 300     N = 400     N = 500\n")
	errortable.write("----------------------------------------------------------------------------\n")
	for i in range(5):
		errortable.write("rho_max = %2d " % X[i])
		for j in range(5):
			errortable.write("|  %.4f   " % Errors[k][i][j])
		errortable.write("\n")
	errortable.write("\n")
errortable.close()

#plot eigenfunctions for one particle in 3D

n = 400
x = 5

A,s = readarrays("Quantum_Eigenpairs_N_%d_X_%d.txt" % (n,x))
X = linspace(0,x,len(A[1]))

plt.figure()
for i in range(1,5):
	print A[0][i-1]
	print sum(abs(A[i])**2)
	plt.plot(X,abs(A[i])**2)

plt.title("Eigenfunctions for one particle \n in 3D harmonic oscillator potential. \n Radial equation, N= %d, $\\rho_{max}$ = %d" % (n,x))
plt.xlabel("$\\rho$")
plt.ylabel("$|u(\\rho)|^2$")
plt.legend(["$\lambda_1 = 3.0000$","$\lambda_2 = 6.9998$","$\lambda_3 = 11.000$","$\lambda_4 = 15.005$"])
#plt.savefig("Eigenfuncs N=%d X=%d" % (n,x))


print "------------------------------"


#plot W-eigenfunctions for one particle in 3D

n = 400
x = [30,10,10,10]
W = [0.01,0.5,1.,5.]


fig,plts = plt.subplots(2,2,sharey="row")
plt.subplots_adjust(hspace=0.3)
fig.suptitle("Groundstate of two electrons in \n 3D harmonic oscillator")
plts = [plts[0][0],plts[0][1],plts[1][0],plts[1][1]]


for i in range(4):
	A,s = readarrays("Quantum_Eigenpairs_1particle_N_%d_X_%d_w_%.2f.txt" % (n,x[i],W[i]))
	X = linspace(0,x[i],len(A[1]))

	#print A[0][0]
	#print sum(abs(A[1])**2)
	plts[i].plot(X,abs(A[1])**2)

	B,t = readarrays("Quantum_Eigenpairs_2particles_N_%d_X_%d_w_%.2f.txt" % (n,x[i],W[i]))
	Xb = linspace(0,x[i],len(B[1]))

	#print B[0][0]
	#print sum(abs(B[1])**2)
	plts[i].plot(Xb,abs(B[1])**2,'r')
	plts[i].legend(['Non-interacting','Interacting'])

	plts[i].set_title("$\omega_r$ = %.2f" % W[i])
	if i >= 2:	
		plts[i].set_xlabel("$\\rho$")
	if i == 0 or i == 2:	
		plts[i].set_ylabel("$|u(\\rho)|^2$")

#plt.savefig("Eigenfuncs_2particles N=%d" % n)


print "----------------------"


#plot N transformations

plt.figure()

Ntrans = readarrays("Transformations and Runtimes.txt")[0]
interpolate = polyfit(Ntrans[0],Ntrans[1],2)
Xtrans = arange(0,2*max(Ntrans[0]))
polynomial = interpolate[0]*Xtrans**2 + interpolate[1]*Xtrans + interpolate[2]

plt.plot(Xtrans,polynomial,Ntrans[0],Ntrans[1],'ro')
plt.title("Interpolated polynomial for \n number of transformations as function of N")
plt.legend(['%.2f x**2 + %.2f x + %.2f' % (interpolate[0],interpolate[1],interpolate[2]),"Calculated values"])
plt.xlabel("N")
plt.ylabel("Number of transformations")
#plt.savefig("Interpol_Transformations.png")

print "--------------------------"

#plot Runtimes

plt.figure()
interpolate = polyfit(Ntrans[0],Ntrans[2],4)
Xtrans = arange(0,2*max(Ntrans[0]))
polynomial = interpolate[0]*Xtrans**4 + interpolate[1]*Xtrans**3 + interpolate[2]*Xtrans**2 + interpolate[3]*Xtrans + interpolate[4]

plt.plot(Xtrans,polynomial,Ntrans[0],Ntrans[2],'ro')
plt.title("Interpolated polynomial for \n Runtimes as function of N")
plt.legend(['%.2f x**4 + %.2f x**3 + %.2f x**2 + %.2f x + %.2f' % (interpolate[0],interpolate[1],interpolate[2],interpolate[3],interpolate[4]),"Calculated values"])
plt.xlabel("N")
plt.ylabel("Runtimes")
print "Runtime Coefficients:" 
print interpolate

#plt.savefig("Interpol_Runtimes.png")

plt.show()
