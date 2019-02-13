# Abbey Johnson
# January 28, 2019
# Math 578, Lab 3: 1D Diffusion - comparison with exact solution
# Python Code

# 1-dimensional diffusion process, with (constant) diffusivity D > 0,
# for the concentration u(x,t) in the semi-infinite interval 0 =< x < infinity
# during some time 0 =< t =< t_end
# initially, concentration is uniform, u_init(x) = 0

# this code is designed for the diffusion equation with Dirichlet boundary conditions

# problem admits explicit similarity solution: u(x,t) = 1 - erf( x / ( 2 sqrt(D t) ) ) = erfc(b / 2 sqrt(D t) ) )

# debugging: test diffusion code (ie, numerical solution) against the exact solution on interval [0, b]

# import numpy for arrays
import numpy as np
# import math module
import math

if __name__ == "__main__":
    print "Lab3 is being run directly \n"
else:
    print "Lab3 is being imported \n"

# open files
init_file = open("/Users/ajohn261/Desktop/Lab3Init.txt", "r")
outf = open('prof.txt', 'w')
comp = open('compare.txt', 'w')

# BEGIN SUBROUTINES

# mesh of M uniform control volumes in the interval [a, b] up to time t_end
def MESH(a, dx, Mp1, b):
    # left boundary of mesh
    x[0] = a
    # internal nodes of mesh
    x[1] = a + dx/2
    for i in range(2, Mp1):
        x[i] = x[1]+(i-1)*dx
    # right boundary of mesh
    x[Mp1] = b
    return x

# initialization subroutine (for the initial concentration profile / initial condition)
def INIT(Mp1):
    for j in range(0, (Mp1+1)):
        U[j] = 0
    return U

# output subroutine
def OUTPUT(time, Mp1, x, U):
    # print x(i), U(i) into output file
    outf.write('# Profile of u(x,t) at time = %f \n' % time)
    for k in range(0, (Mp1+1)):
        outf.write('%f %f \n' % (x[k], U[k]))
    outf.write('\n')

# flux subroutine
def FLUX(Mp1, D, U, x, b, time):
    # left boundary condition (Dirichlet)
    U[0] = 1.
    # right boundary condition (Dirichlet) (caution: time-dependent boundary condition)
    # right boundary condition, impose exact solution at x = b
    U[Mp1] = math.erfc(b/(2*math.sqrt(D*time)))
    # left boundary flux
    # F[1] = -(D*(U[1]-U[0])/(x[1]-x[0]))
    # flux at internal faces
    for ii in range(1, (Mp1+1)):
        F[ii] = -D * (U[ii]-U[ii-1]) / (x[ii]-x[ii-1])
    # right boundary flux
    # F[Mp1] = -(D*(U[Mp1]-U[M]))/(x[Mp1]-x[M])
    return F

# PDE subroutine
def PDE(Mp1, dt, dx, F, U, b, D, time):
    # internal concentration
    for jj in range(1, Mp1):
        U[jj] = U[jj]+(dt/dx)*(F[jj]-F[jj+1])
    # update boundary concentrations after internal concentrations have been updated
    # left boundary concentration
    # U[0] = U[1]+(F[1]/D)*(x[1]-x[0])
    # left Dirichlet condition
    U[0] = 1.
    # right boundary concentration
    # U[Mp1] = U[M]-(F[Mp1]/D)*(x[Mp1]-x[M])
    # right Dirichlet condition (caution: time-dependent) (exact solution at x = b)
    U[Mp1] = math.erfc(b/(2*math.sqrt(D*time)))
    return U

# Trapezoidal Rule Quadrature
def TrapzRule(Mp1, dx, U):
    Area = (U[0]-U[1]-U[M]+U[Mp1])/4
    for kk in range(1, Mp1):
        Area = Area + U[kk]
    quad = Area * dx
    return quad

def COMPARISON(x, U, Mp1, time, D):
    comp.write('# x(i) U(i) u_exact(i) error(i) at time = %f \n' % time)
    ERR = 0
    for q in range(0, (Mp1+1)):
        # compute exact solution
        u_exact[q] = math.erfc(x[q]/(2*math.sqrt(D*time)))
        # compute error at x[i]
        ERRi = abs(U[q]-u_exact[q])
        # print x(i), U(i) into comparison output file
        comp.write('%f %f %f %f \n' % (x[q], U[q], u_exact[q], ERRi))
        # compute max error
        ERR = max(ERRi, ERR)
    comp.write('\n')
    return ERR

# END SUBROUTINES

# BEGIN TIME-STEPPING SCHEME

# read from data file: MM, tend, factor, dtout, D, a, b
contents = init_file.read().split()
MM = np.float64(contents[0])
tend = np.float64(contents[1])
factor = np.float64(contents[2])
dtout = np.float64(contents[3])
D = np.float64(contents[4])
a = np.float64(contents[5])
b = np.float64(contents[6])

dx = np.float64(1.0/MM)

# M = int((b-a)/dx)
M = int((b-a)*MM)

Mp1 = M+1

dtEXPL = np.float64((dx*dx)/(2*D))

# time-step (for stability of the explicit scheme)
dt = np.float64(factor*dtEXPL)

# declare variables and dimensions of arrays
x = np.zeros((M+2), dtype=np.float64)
U = np.zeros((M+2), dtype=np.float64)
F = np.zeros((M+2), dtype=np.float64)
u_exact = np.zeros((M+2), dtype=np.float64)

# call mesh subroutine
x = MESH(a, dx, Mp1, b)

# initialize
nsteps = 0
time = 0.0
tout = dtout
MaxSteps = int(tend/dt)+1

# call initialization subroutine
U = INIT(Mp1)

# call output routine to record the concentration profile at time t = 0
# OUTPUT(time, Mp1, x, U)

quad = TrapzRule(Mp1, dx, U)
print('Trapezoidal Rule Quadrature at time %f is %f \n' % (time, quad))

# begin time-stepping
for nsteps in range(1, (MaxSteps+1)):
    # update time = time + dt
    time = nsteps * dt
    # call flux subroutine
    F = FLUX(Mp1, D, U, x, b, time)
    # call PDE subroutine
    U = PDE(Mp1, dt, dx, F, U, b, D, time)
    if time >= tout:
        # call compare subroutine (as long as exact solution is available)
        ERR = COMPARISON(x, U, Mp1, time, D)
        # call output subroutine
        OUTPUT(time, Mp1, x, U)
        # update tout
        tout = tout + dtout
    # end if-statement
# end time-stepping for-loop

# print run-time information to screen
print('DONE: exiting at t = %f after %i steps. \n' % (time, nsteps))

quad = TrapzRule(Mp1, dx, U)
print('Trapezoidal Rule Quadrature at time %f is %f \n' % (time, quad))

# print maximum error at t_end for i=0:M+1
print('Maximum error is %e \n' % ERR)

print('Parameters for this run: ')
print('MM = %f ' % MM)
print('t_end = %f ' % tend)
print('factor = %f ' % factor)
print('dtout = %f ' % dtout)
print('D = %f ' % D)
print('a = %f ' % a)
print('b = %f \n' % b)

# close files
init_file.close()
outf.close()
comp.close()

# END TIME-STEPPING SCHEME
