using LinearAlgebra
export G_0
export G_0_near

#x polarized
export G_0_xx
export G_0_near_xx


#------------------------------------------------------------------------------#
#define the free space Green's function, divergence at the origin ignored
function G_0(r,k)
	rho=norm(r)
	A=exp(im*k*rho)/(4*pi*(k^2)*(rho^3)) 

	c1=(k*rho)^2 + im*k*rho - 1.0
	c2=-(k*rho)^2 - 3.0*im*k*rho + 3.0

	return A*(  c1*I + c2*( r*transpose(r) )/(rho^2)   )
end
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#define the near field Green's function, divergence at the origin ignored
function G_0_near(r,k)
	rho=norm(r)
	A=1.0/(4*pi*(k^2)*(rho^3)) 
	return A*(  -I + (3.0/rho^2)*( r*transpose(r) )   )
end
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#define the free space Green's function projected along x, divergence at the origin ignored
function G_0_xx(r,k)
	rho=norm(r)
	A=exp(im*k*rho)/(4*pi*(k^2)*(rho^3)) 

	c1=(k*rho)^2 + im*k*rho - 1.0
	c2=-(k*rho)^2 - 3.0*im*k*rho + 3.0

	return A*(  c1 + c2*(r[1]/rho)^2   )
end
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#define the free space Green's function projected along x, divergence at the origin ignored
function G_0_near_xx(r,k0)
	rho=norm(r)
	A=1.0/(4*pi*(k0^2)*(rho^3)) 
	return A*(  -1.0 + 3.0*(r[1]/rho)^2   )
end
#------------------------------------------------------------------------------#
