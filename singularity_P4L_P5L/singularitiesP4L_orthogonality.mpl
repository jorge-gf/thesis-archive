# ==============================================================================
# SINGULARITY ANALYSIS FOR P4L with ORTHOGONALITY & PARALLELISM constraints
# files for PhD thesis Singularity and Stability Analysis of vision-based controllers
# JORGE GARCIA FONTAN
# ==============================================================================
restart;

# LOAD packages and functions
with(LinearAlgebra): with(Groebner): with(CodeGeneration):

# Function: Compute the square submatrices of size p of A
Minors_submatrices:=[(A::Matrix,p::posint)->
seq(seq(A(r,c), c=combinat:-choose([$1..op(1,A)[2]],p)),
r=combinat:-choose([$1..op(1,A)[1]],p))]:

# Function: Determinant minors of order p of A
Minors_determinants:=[(A::Matrix,p::posint)->
seq(seq(LinearAlgebra:-Determinant(A(r,c)),
c=combinat:-choose([$1..op(1,A)[2]],p)),
r=combinat:-choose([$1..op(1,A)[1]],p))]:

# ==============================================================================

# PARAMETRIZATION:
# A line is described by a point and a direction, used to compute the Plucker coordinates

OM1 := <0,0,0>:     CM1 := -<x,y,z> + OM1:
OM2 := <0,0,d1>:    CM2 := -<x,y,z> + OM2:
OM3 := <d2,d3,0>:   CM3 := -<x,y,z> + OM3:
OM4 := <0,d4,d5>:   CM4 := -<x,y,z> + OM4:

# Direction of lines
U1 := <1,0,0>:        L1 := (U1 &x CM1):  # Plucker coordinates
U2 := <0,1,0>:        L2 := (U2 &x CM2):
U3 := <0,0,1>:        L3 := (U3 &x CM3):
U4 := <1,0,0>:        L4 := (U4 &x CM4):

Ux1 := U1[1]: Uy1 := U1[2]: Uz1 := U1[3]:
Ux2 := U2[1]: Uy2 := U2[2]: Uz2 := U2[3]:
Ux3 := U3[1]: Uy3 := U3[2]: Uz3 := U3[3]:
Ux4 := U4[1]: Uy4 := U4[2]: Uz4 := U4[3]:

Lx1 := L1[1]: Ly1 := L1[2]: Lz1 := L1[3]:
Lx2 := L2[1]: Ly2 := L2[2]: Lz2 := L2[3]:
Lx3 := L3[1]: Ly3 := L3[2]: Lz3 := L3[3]:
Lx4 := L4[1]: Ly4 := L4[2]: Lz4 := L4[3]:

DELTA1 := sqrt(Lx1*Lx1 + Ly1*Ly1):
DELTA2 := sqrt(Lx2*Lx2 + Ly2*Ly2):
DELTA3 := sqrt(Lx3*Lx3 + Ly3*Ly3):
DELTA4 := sqrt(Lx4*Lx4 + Ly4*Ly4):

# ==============================================================================
(*
# INTERACTION MATRIX (not necessary to compute)
M11 := Transpose(< (Lx1*Ly1*Uz1), (Ly1*Ly1*Uz1), -(Ly1*( Lx1*Ux1 + Ly1*Uy1 )), Lx1*Ly1*Lz1, Ly1*Ly1*Lz1, -Ly1*DELTA1*DELTA1 >)*(1/(DELTA1*DELTA1*DELTA1)):
M12 := Transpose(< -(Lx1*Lx1*Uz1), -(Lx1*Ly1*Uz1), Lx1*( Lx1*Ux1 + Ly1*Uy1 ), -Lx1*Lx1*Lz1, -Lx1*Ly1*Lz1, Lx1*DELTA1*DELTA1 >)*(1/(DELTA1*DELTA1*DELTA1)):
M13 := Transpose(< (Uy1*DELTA1*DELTA1 + Ly1*Lz1*Uz1), -( Ux1*DELTA1*DELTA1 + Lx1*Lz1*Uz1 ), Lz1*( Lx1*Uy1 - Ly1*Ux1 ), Ly1*(Lz1*Lz1 + DELTA1*DELTA1), -Lx1*(Lz1*Lz1 + DELTA1*DELTA1), 0 >)*(1/(DELTA1*DELTA1*DELTA1)):
M1 := <M11; M12 ; M13>:

M21 := Transpose(< (Lx2*Ly2*Uz2), (Ly2*Ly2*Uz2), -(Ly2*( Lx2*Ux2 + Ly2*Uy2 )), Lx2*Ly2*Lz2, Ly2*Ly2*Lz2, -Ly2*DELTA2*DELTA2 >)*(1/(DELTA2*DELTA2*DELTA2)):
M22 := Transpose(< -(Lx2*Lx2*Uz2), -(Lx2*Ly2*Uz2), Lx2*( Lx2*Ux2 + Ly2*Uy2 ), -Lx2*Lx2*Lz2, -Lx2*Ly2*Lz2, Lx2*DELTA2*DELTA2 >)*(1/(DELTA2*DELTA2*DELTA2)):
M23 := Transpose(< (Uy2*DELTA2*DELTA2+ Ly2*Lz2*Uz2), -( Ux2*DELTA2*DELTA2 + Lx2*Lz2*Uz2 ), Lz2*( Lx2*Uy2 - Ly2*Ux2 ), Ly2*(Lz2*Lz2 + DELTA2*DELTA2), -Lx2*(Lz2*Lz2 + DELTA2*DELTA2), 0 >)*(1/(DELTA2*DELTA2*DELTA2)):
M2 := <M21; M22 ; M23>:

M31 := Transpose(< (Lx3*Ly3*Uz3), (Ly3*Ly3*Uz3), -(Ly3*( Lx3*Ux3 + Ly3*Uy3 )), Lx3*Ly3*Lz3, Ly3*Ly3*Lz3, -Ly3*DELTA3*DELTA3 >)*(1/(DELTA3*DELTA3*DELTA3)):
M32 := Transpose(< -(Lx3*Lx3*Uz3), -(Lx3*Ly3*Uz3), Lx3*( Lx3*Ux3 + Ly3*Uy3 ), -Lx3*Lx3*Lz3, -Lx3*Ly3*Lz3, Lx3*DELTA3*DELTA3 >)*(1/(DELTA3*DELTA3*DELTA3)):
M33 := Transpose(< (Uy3*DELTA3*DELTA3+ Ly3*Lz3*Uz3), -( Ux3*DELTA3*DELTA3 + Lx3*Lz3*Uz3 ), Lz3*( Lx3*Uy3 - Ly3*Ux3 ), Ly3*(Lz3*Lz3 + DELTA3*DELTA3), -Lx3*(Lz3*Lz3 + DELTA3*DELTA3), 0 >)*(1/(DELTA3*DELTA3*DELTA3)):
M3 := <M31; M32 ; M33>:

M41 := Transpose(< (Lx4*Ly4*Uz4), (Ly4*Ly4*Uz4), -(Ly4*( Lx4*Ux4 + Ly4*Uy4 )), Lx4*Ly4*Lz4, Ly4*Ly4*Lz4, -Ly4*DELTA4*DELTA4 >)*(1/(DELTA4*DELTA4*DELTA4)):
M42 := Transpose(< -(Lx4*Lx4*Uz4), -(Lx4*Ly4*Uz4), Lx4*( Lx4*Ux4 + Ly4*Uy4 ), -Lx4*Lx4*Lz4, -Lx4*Ly4*Lz4, Lx4*DELTA4*DELTA4 >)*(1/(DELTA4*DELTA4*DELTA4)):
M43 := Transpose(< (Uy4*DELTA4*DELTA4+ Ly4*Lz4*Uz4), -( Ux4*DELTA4*DELTA4 + Lx4*Lz4*Uz4 ), Lz4*( Lx4*Uy4 - Ly4*Ux4 ), Ly4*(Lz4*Lz4 + DELTA4*DELTA4), -Lx4*(Lz4*Lz4 + DELTA4*DELTA4), 0 >)*(1/(DELTA4*DELTA4*DELTA4)):
M4 := <M41; M42 ; M43>:

M_P4L := < M1; M2; M3; M4 >: # Full interaction matrix
*)

# ==============================================================================

# BASIS for the rows of the interaction matrix M_P4L
f11 := simplify(U1 &x CM1):     m12 := simplify(U1 &x f11):
f21 := simplify(U2 &x CM2):     m22 := simplify(U2 &x f21):
f31 := simplify(U3 &x CM3):     m32 := simplify(U3 &x f31):
f41 := simplify(U4 &x CM4):     m42 := simplify(U4 &x f41):

xi11 := Vector([f11, OM1 &x f11]): xi12 := Vector([0,0,0,m12]):
xi21 := Vector([f21, OM2 &x f21]): xi22 := Vector([0,0,0,m22]):
xi31 := Vector([f31, OM3 &x f31]): xi32 := Vector([0,0,0,m32]):
xi41 := Vector([f41, OM4 &x f41]): xi42 := Vector([0,0,0,m42]):
XI := Matrix([xi11, xi21, xi31, xi41, xi12, xi22, xi32, xi42]):

XI_minors := Minors_determinants(XI, 6): # Maximal minors of the basis XI
XI_minors := factor(convert(convert(XI_minors, set)  minus {0}, list)): # eliminates redundant and null elements
XI_minors := expand(XI_minors):

# ==============================================================================

# CONDITIONS FOR DEGENERACY
# hyperbolic equations
force_triplet := combinat:-choose(map(evala, [f11, f21, f31, f41]), 3):
Gijk := []:
for ii from 1 to numelems(force_triplet) do
  Gijk := [op(Gijk), Transpose(force_triplet[ii][1]).( force_triplet[ii][2] &x force_triplet[ii][3] )]:
end do:
Gijk := collect(Gijk, [x,y,z], distributed):

# cubic equations - no longer all cubic for orthogonal lines
moment_triplet := combinat:-choose([m12, m22, m32, m42], 3):
Hlmn := []:
for ii from 1 to numelems(moment_triplet) do
  Hlmn := [op(Hlmn), Transpose(moment_triplet[ii][1]).( moment_triplet[ii][2] &x moment_triplet[ii][3] )]:
end do:
Hlmn := collect(Hlmn, [x,y,z], distributed):

I16 := []:     # Ideal formed by polynomials Gijk*Hlmn
for i from 1 to 4 do
  for j from 1 to 4 do
    I16 := [op(I16), Gijk[i]*Hlmn[j]];
  end do:
end do:

XI_minors := expand(XI_minors):   # All the minors of the basis
I16 := expand(I16):

# Remaining 6 minors of degree 5 - obtained by removing the terms Gijk*Hlmn from the full set
K6 := convert(convert(XI_minors,set) minus convert(I16, set) minus convert(-I16, set), list): # Extract the remaining 6 minors

# The singularities are given by V(XI_minors) = V(Gijk, K6) U V(Hlmn, K6)
# i.e. the common roots of XI_minors


# ==============================================================================
# EVALUATE VARIETY V(Gijk, K6) - Common roots of (Gijk & K6)

# Compute Grobner Basis of Gijk wrt an ordering grevlex(z,y,x, d1,d2,d3,d4,d5)
# since the parameters are considered as variables, the result is a GB under any specialization of the parameters
gb_G_full := Groebner:-Basis(Gijk, tdeg(z,y,x, d1,d2,d3,d4,d5)):
gb_G := Groebner:-Basis(Gijk, plex(x,y,z)):

# Compute Normal Form of the polynomials in K6 by the elements of the GB
NF_K6_G := []:  # Normal Form of the polynomials (i.e. remainder of polynomial division)
quot_NF_K6_G := []:   # Quotients of the polynomial division
for i from 1 to numelems(K6) do
  NF := Groebner:-NormalForm(K6[i], gb_G_full, tdeg(z,y,x, d1,d2,d3,d4,d5), quot):
  NF_K6_G := [op(NF_K6_G), NF]:
  quot_NF_K6_G := [op(quot_NF_K6_G), quot]:
NF := 'NF':
quot := 'quot':
od:

# Note that NF_K6_G contains only zeroes - therefore V(Gijk, K6) = V(Gijk),
# and therefore V(Gijk) is a component of the singularity loci

factor(gb_G[1]);  # gb_G[1] factorises in two linear terms - defining two planes

# by substituting z=d1 and d5*y=d4*z in the elements of gb_G we obtain the equations for the two transversal lines

transversal1 := Groebner:-Basis([op(Gijk), -z + d1], plex(x,y,z)):
transversal2 := Groebner:-Basis([op(Gijk), d4*z - d5*y], plex(x,y,z)):

# ==============================================================================
# EVALUATE VARIETY V(Hlmn, K6) - Common roots of (Hlmn & K6)

# Compute Grobner Basis of Hlmn wrt an ordering grevlex(z,y,x, d1,d2,d3,d4,d5)
gb_H := Groebner:-Basis(expand(Hlmn), tdeg(z,y,x)):
gb_H_full := Groebner:-Basis(Hlmn, tdeg(z,y,x, d1,d2,d3,d4,d5)):

# Compute Normal Form of the polynomials in K6 by the elements of the GB
quot_NF_K6_H := []:
NF_K6_H := []:
for i from 1 to numelems(K6) do
  NF := Groebner:-NormalForm(K6[i], gb_H_full, tdeg(z,y,x, d1,d2,d3,d4,d5), quot):
  NF_K6_H := [op(NF_K6_H), NF]:
  quot_NF_K6_H := [op(quot_NF_K6_H), quot]:
NF := 'NF':
quot := 'quot':
od:

# The terms NF_K6_H are not zero, hence we have V(Hlmn, K6) \neq V(Hlmn)

# GB of (Hlmn,K6)
gb_H_K := Groebner:-Basis([op(Hlmn), op(K6)], tdeg(z,y,x, d1,d2,d3,d4,d5)):

# There are 4 solutions to the system (Hlmn,K6). We verify that they all lie on the observed lines
# hence there are no isolated point singularities
sols_H_K := solve(gb_H_K, [x,y,z]):

# ==============================================================================
# END
# ==============================================================================
