#routines for Spin Model numerics in the single excitation manifold
using LinearAlgebra
using StaticArrays

#include the Green's function definition
include("./GreensTensor.jl")


#------------------------------------------------------------------------------#
#build the Hamiltonian at 1st excitation for the Spin Model
#for polarized atoms along P
#it needs:
# z: an array built by the lattice function
# k₀: resonance momentum of the atom
# Pa: the polarization vector of atoms
function H_1e(z,k₀,P)
 N=length(z)
 H=zeros(Complex{Float64},N,N)
 for m=1:N,l=1:N
    if l!=m
      H[l,m] =  -(3pi/k₀)*dot(P,G_0(z[l]-z[m],k₀)*P)
    else end
 end
 return H
end
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#Electrical field scattered from atoms
# - in a point r (3D vector)
# - of a given atomic state s (column vector!!)
# - of a lattice z (lattice object)
# - polarization of the atoms P
# - transition frequency of the atom k₀
function E_scatt(r,z,P,c,k₀)
  E_plus=zeros(Complex{Float64},3)
  for l in 1:length(z)
    E_plus .+= (3pi/k₀)*c[l]*(G_0(r-z[l],k₀)*P)
  end
  return E_plus
end
#------------------------------------------------------------------------------#
