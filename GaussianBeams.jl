#Module for creating Various Gaussian beams
using LinearAlgebra

import SpecialFunctions.besselj0
import SpecialFunctions.besselj1
import SpecialFunctions.erfi

using QuadGK

#------------------------------------------------------------------------------#
#Gaussian beam in the paraxial approx propagating along z
#    r: point at which we evaluate the field
#  E_0: amplplitude
#    p: polarization transverse x,y) 3d vector
#  w_0: beam waist
#  k: 3d wavevector, typically with norm the resonant k0
function E_Laguerre_Gauss(r,E_0,p,vk,w_0)

  #auxiliary variables for the calculation
  k=norm(vk)
  zk=dot(vk,r)/k
  rho2=norm(r)^2-zk^2

  z_R=0.5*k*w_0^2
  w=w_0*sqrt(1.0+(zk/z_R)^2)
  R=zk*(1+(z_R/zk)^2)
  phi=atan(zk/z_R)
  
  #if zk==0 the calculation simplifies 
  if zk==0
    return p*E_0*exp(-rho2/(w_0^2))
  else
    return (p*E_0*w_0/w)*exp( im*(k*zk + 0.5*k*rho2/R - phi) )*exp(-rho2/(w^2))
  end
end
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
function norm_E_Laguerre_Gauss(E0,w0)
  #norm of the field integrated over a 2D plane
  return 0.5*pi*(E0*w0)^2
end
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#Gaussian beam (exact solution)  propagating along z
#    r: point at which we evaluate the field
#  E_0: amplplitude
#    p: polarization 3d vector
#  w_0: beam waist
#  k: 3d wavevector, typically with norm the resonant k0
#it is not enforced that k and p are orthogonal
function E_Gauss(r,E_0,p,vk,w_0)

  zk=dot(vk,r)
  rho=sqrt(norm(r)^2-zk^2)
  k=norm(vk)

  f_p(b)=b*exp( -(b*k*w_0/2.0)^2 )*exp(im*zk*sqrt(1.0-b^2))*besselj0(b*k*rho)
  res_p=quadgk(f_p,0.0,1.0;atol=1.0e-10)

#  #along the polarization vector
  rp=dot(p,r)
  res_k=[im*0.0,im*0.0]
  if rp!=0.0 && rho!=0
   #change of variable respect to Marco's paper b=sin(s)    
   f_k(s)=(rp/rho)*exp(im*zk*cos(s))*(sin(s)^2)*exp( -(sin(s)*k*w_0/2.0)^2 )*besselj1(sin(s)*k*rho)
   res_k=quadgk(f_k,0,pi/2.0;atol=1.0e-10)
  else end

#  #polarization component along the propagation
  pk=vk/k  
  return E_0*( (res_p[1]*p) - (im*res_k[1]*pk))
end

function norm_E_Gauss(E_0,k0,w_0)
  #norm of the field integrated over a 2D plane
  λ0=2*pi/k0
  beta=k0*w_0

  #if the beam wait is sufficiently large the paraxial approx applies
  if (w_0/λ0)>4.5
    return 2*pi/(beta^2)
  else
    D=0.5*sqrt(pi)*exp(-0.5*beta^2)*erfi(beta/sqrt(2.0))
    return (pi*(E_0^2)/(beta^2))*(1.0+sqrt(2.0)*(-beta^-1 + beta)*D)
  end
end
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
#=
Gaussian field at atomic positions
z: list of atomic positions (3D vectors)
Pa: polarization of the atoms
E_0: field amplitude
Pf: polarization of the field
vk: field wavevector
w_0: beam waist
fieldtype: string to specify the field type
=#
function E_at_atoms(z,Pa,E_0,Pf,vk,w_0,fieldtype)

  if fieldtype=="Gauss"
    return [dot(Pa,E_Gauss(p,E_0,Pf,vk,w_0)) for p in z]       
  elseif fieldtype=="Laguerre-Gauss"
    return [dot(Pa,E_Laguerre_Gauss(p,E_0,Pf,vk,w_0)) for p in z]      
  else printstyled("Error: $fieldtype is not a valid field type.\n"; color=:red)
   end
   
end
#------------------------------------------------------------------------------#

