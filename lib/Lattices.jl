#=
Module to create arrays of atomic positions in 3D, with the desired Bravais lattice.
There are 14 possible configurations. Here we cover
Crystal family: Cubic
    TODO: Triclinic,Monoclinic,Orthorhombic,Tetragonal,Hexagonal,
Centering type: Primitive (P), Body-centered (I), Face-centered (F)
    TODO: Base-centered (B)

Arguments:
  l0: center of the lattice
  d: lattice constant
  N_x,N_y,N_z: number of positions to return in each dimension
Return:
  L: Array of 3D StaticArrays
=#

using LinearAlgebra
using StaticArrays

#------------------------------------------------------------------------------#
function CubicP(l0,d,N_x,N_y,N_z)
  #starting point
  z0=SVector{3,Float64}(l0[1]-0.5*d*(N_x-1),l0[2]-0.5*d*(N_y-1),l0[3]-0.5*d*(N_z-1))
  #allocate
  L=Array{SVector{3,Float64},1}(undef,0)
  for i=1:N_x, j=1:N_y, k=1:N_z
    step=d*SVector{3,Float64}(i-1,j-1,k-1)   
	  push!(L,z0+step)
  end
  return L
end
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
function CubicI(l0,d,N_x,N_y,N_z)
   #starting point
   z0=SVector{3,Float64}(l0[1]-0.5*d*(N_x-1),l0[2]-0.5*d*(N_y-1),l0[3]-0.5*d*(N_z-1))
   #displacement
   v=d*SVector{3,Float64}(0.5,0.5,0.5)
   return vcat(CubicP(l0,d,N_x,N_y,N_z),CubicP(l0+v,d,N_x,N_y,N_z))
end
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
function CubicF(l0,d,N_x,N_y,N_z)
   #starting point
   z0=SVector{3,Float64}(l0[1]-0.5*d*(N_x-1),l0[2]-0.5*d*(N_y-1),l0[3]-0.5*d*(N_z-1))
   #displacements
   vxy=d*SVector{3,Float64}(0.5,0.5,0.0)
   vxz=d*SVector{3,Float64}(0.5,0.0,0.5)
   vyz=d*SVector{3,Float64}(0.0,0.5,0.5)

   return vcat(CubicP(l0,d,N_x,N_y,N_z),CubicP(l0+vxy,d,N_x,N_y,N_z),CubicP(l0+vxz,d,N_x,N_y,N_z),CubicP(l0+vyz,d,N_x,N_y,N_z))
end
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
function Diamond(l0,d,N_x,N_y,N_z)
   D_half=SVector{3,Float64}(d/8,d/8,d/8)
   return vcat(CubicF(l0-D_half,d,N_x,N_y,N_z),CubicF(l0+D_half,d,N_x,N_y,N_z))
end
#------------------------------------------------------------------------------#
