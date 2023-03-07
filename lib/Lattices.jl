using LinearAlgebra
using StaticArrays

#------------------------------------------------------------------------------#
function lattice_for_plot(l)
 N=length(l)
 
 X = [l[i][1] for i=1:N]
 Y = [l[i][2] for i=1:N]
 Z = [l[i][3] for i=1:N]
 
 return [X Y Z]
end
#------------------------------------------------------------------------------#


#==============================================================================#

#-------------------14 Bravais lattices----------------------------------------#

#Crystal family:Triclinic,Monoclinic,Orthorhombic,Tetragonal,Hexagonal,Cubic	
#Primitive (P) 	Base-centered (B) 	Body-centered (I) 	Face-centered (F)

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


#------------------------------------------------------------------------------#
#-----------------------RECIPROCAL LATTICES------------------------------------#
#------------------------------------------------------------------------------#
reciprocal_CubicP(d,N_x,N_y,N_z)=CubicP([0.0,0.0,0.0],2pi/d,2N_x-1,2N_y-1,2N_z-1)
reciprocal_CubicI(d,N_x,N_y,N_z)=CubicF([0.0,0.0,0.0],4pi/d,2N_x-1,2N_y-1,2N_z-1)
reciprocal_CubicF(d,N_x,N_y,N_z)=CubicI([0.0,0.0,0.0],4pi/d,2N_x-1,2N_y-1,2N_z-1)
#------------------------------------------------------------------------------#




#------------------------------------------------------------------------------#
#---------------------PATHS IN BRILLOUIN ZONE---------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#the function create a path between the points a and b; this 3d path can be
#followed accessing every 3d vector on the second index v_i=v[:,i], the first
#index are the coordinates x,y,z
function path_3d(a,b,N)
 x=range(a[1],stop=b[1],length=N)
 y=range(a[2],stop=b[2],length=N)
 z=range(a[3],stop=b[3],length=N)
 return [SVector{3,Float64}(x[i],y[i],z[i]) for i=1:N]
end
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#----------IRREDUCIBLE PATHS for square/cubic lattices-------------------------#
#The following functions build an irreducible path of momentums in the first
#Brilloin zone, to explore the band structure of a square/cubic lattice
#------------------------------------------------------------------------------#
##These functions just needs the lattice constant and the number of steps that the
##user whants to create to cover an irreducible path

#------------------------------------------------------------------------------#
#-----Irreducible path in 2D---------------------------------------------------#

function irr_path_Square_XY(d,n)

  X=[pi/d,0.0,0.0]
  G=[0.0,0.0,0.0]
  M=[pi/d,pi/d,0.0]
  path1=path_3d(X,G,n)
   deleteat!(path1,n)
  path2=path_3d(G,M,n)
   deleteat!(path2,n)
  path3=path_3d(M,X,n)
  return vcat(path1,path2,path3)

end
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#-----Irreducible path in 3D---------------------------------------------------#
function irr_path_CubicP(d,n)
  Γ=[0.0,0.0,0.0]
  X=[pi/d,0.0,0.0]
  M=[pi/d,pi/d,0.0]
  R=[pi/d,pi/d,pi/d]
  path1=path_3d(Γ,X,n)
   deleteat!(path1,n)
  path2=path_3d(X,M,n)
   deleteat!(path2,n)
  path3=path_3d(M,R,n)
   deleteat!(path3,n)
  path4=path_3d(R,Γ,n)
  return vcat(path1,path2,path3,path4)
end
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
function irr_path_CubicF(d,n)

  X=[0,2π/d,0]
  U=[π/2d,2π/d,π/2d]
  L=[π/d,π/d,π/d]
  Γ=[0,0,0]
  W=[π/d,2π/d,0]
  K=[3π/2d,3π/2d,0] 

  path1=path_3d(X,U,n)
   deleteat!(path1,n)
  path2=path_3d(U,L,n)
   deleteat!(path2,n)
  path3=path_3d(L,Γ,n)
   deleteat!(path3,n)
  path4=path_3d(Γ,X,n)
   deleteat!(path4,n)
  path5=path_3d(X,W,n)
   deleteat!(path5,n)
  path6=path_3d(W,K,n)
  
  return vcat(path1,path2,path3,path4,path5,path6)
end
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
function irr_path_CubicI(d,n)

  Γ=[0,0,0]
  H=[0,0,2π/d]  
  P=[π/d,π/d,π/d]  
  N=[0,π/d,π/d]  

  path1=path_3d(Γ,H,n)
  path2=path_3d(H,P,n)
  path3=path_3d(P,N,n)
  path4=path_3d(N,Γ,n)
  return vcat(path1,path2,path3,path4)
end
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#----------------------LATTICES WITH SHAPE-------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
function CubicP_sphere(d,R)
  #number of atoms within radius
  N=floor(Int,R/d)
   
  #allocate
  L=Array{SVector{3,Float64},1}(undef,0)
  for i=-N:N, j=-N:N, k=-N:N
	new_point=d*SVector{3,Float64}(i,j,k)
        if norm(new_point)<R   
	  push!(L,new_point)
	else end
  end
  return L
end
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
function CubicP_cylinder(d,R,Nz)
  #number of atoms within radius
  N=floor(Int,R/d)
  
  z0=SVector{3,Float64}(0.0,0.0,-0.5*d*(Nz-1))
  
  #allocate
  L=Array{SVector{3,Float64},1}(undef,0)
  for i=-N:N, j=-N:N, k=1:Nz
	new_point=d*SVector{3,Float64}(i,j,k-1)
	norm_xy=norm(new_point[1:2])
        if norm_xy < 1.01*R   
	  push!(L,z0+new_point)
	else end
  end
  return L
end
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#----------------------LATTICES WITH HOLES-------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
function One_Punch_Square_Hole_xy(lat,Lh)
  #initialize a new empty array
  L=Array{SVector{3,Float64},1}(undef,0)

  #number of elements in the array
  n=length(lat)
  for i=1:n
    p=lat[i]
    c1=abs(p[1])>1.01*Lh
    c2=abs(p[2])>1.01*Lh
    if c1 || c2
      push!(L,p)
    else end
  end
  return L
end
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#extract the indeces of the atoms at the border of the hole
function ind_inner_border(lat,Lh)
 ind=Array{Int,1}(undef,0)
 for i=1:length(lat)
  cond1=abs(lat[i][1])<(Lh+1.01)
  cond2=abs(lat[i][2])<(Lh+1.01)
  if cond1 && cond2
   push!(ind,i)
  else end
 end
 return ind
end
#------------------------------------------------------------------------------#
