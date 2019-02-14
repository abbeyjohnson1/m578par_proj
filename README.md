# Math 578: Parallelization Project

Copied from Lab 3:

Consider a (1-dimensional) diffusion process, with (constant) diffusivity D > 0, for the concentration u(x,t) in the semi-infinite interval 0 ≤ x < ∞, during some time  0 ≤ t ≤ tend. Initialy, the concentration is uniform: uinit(x)=0, and the endpoint at x=0 is raised and held at concentration u0(t)=1. 

This problem admits the explicit similarity solution:

		u(x,t) = 1 − ERF( x / (2 SQRT(D t)) ) ,

with ERF(z) the error function.  

To test your diffusion code against this exact solution on an interval [0,b], 
need to impose the exact solution at x=b:

	U(M+1) = 1.0−ERF( b / (2*SQRT(D*t)) ) == ERFC( b / (2*SQRT(D*t)) ) .

So here both BCs are Dirichlet type, and the one at x=b is time-dependent...
think when time needs to be updated...
