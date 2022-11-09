#=================================================================================
# COMPUTATION OF CRITICAL POINTS OF IBVS FROM 4 POINTS
#=================================================================================

restart;
with(LinearAlgebra): with(CodeGeneration): with(Groebner):
interface(rtablesize=16);

Round := (x,n)-> parse(sprintf(sprintf("%%.%df",n),x));

k := 4; # NUMBER OF POINTS
#=================================================================================

# LOAD PARAMETERS for different configurations
# Parameters are :
# d = [d12,d13,....]: vector of distances (squared between each two points)
# sd = [x1_star, ..., y1_star]: vector of desired visual features
# file_param contains also info such as:
# point coordinates in F0-frame
# desired camera position : t_star (position vector)
# desired camera orientation : theta_star, tu_star (tu-vector)

exp_name := 'genericPlanar1':  # select experiment name/number
file_param := cat("/dsk/l1/jorgegf/systemsThesis/parameters/", exp_name, ".mpl"):
read(file_param);   # Import system parameters

Rot_star := Student[LinearAlgebra]:-RotationMatrix(theta_star, tu_star); # desired camera rotation matrix

# OUTPUT FILES: msolve input, output and log files
fname1 := cat("/dsk/l1/jorgegf/4points/october/systems/", exp_name, ".original.sys"):
fname2 := cat("/dsk/l1/jorgegf/4points/october/sols/", exp_name, ".original.sols"):
fname_log := cat("/dsk/l1/jorgegf/4points/october/logs/", exp_name, ".original.log"):

(*
# --------------------------------------------------------------------------------
# (Alternatively) DEFINE PARAMETERS MANUALLY and name files
d := []: # Vector of distances - 6 coordinates
sd := []: # Vector of desired features - 8 coordinates

fname1 := "system.in":
fname2 := "system.out":
fname_log := "system.log":
# --------------------------------------------------------------------------------
*)

d12 := d[1]: d13 := d[2]: d14 := d[3]: d23 := d[4]: d24 := d[5]: d34 := d[6]:

xi_star := Vector(4): yi_star := Vector(4):
xi_star[1] := sd[1]: xi_star[2] := sd[2]: xi_star[3] := sd[3]: xi_star[4] := sd[4]: yi_star[1] := sd[5]: yi_star[2] := sd[6]: yi_star[3] := sd[7]: yi_star[4] := sd[8]:

#=================================================================================
# VARIABLES IN (s,Z)-space
xi := <x1,x2,x3,x4>: yi := <y1,y2,y3,y4>:
Zi := <Z1,Z2,Z3,Z4>:
#s := [x1,y1,x2,y2,x3,y3,x4,y4]:

Xi := xi*~Zi: Yi := yi*~Zi:

m1 := <Xi[1],Yi[1],Zi[1]>: # # POINT COORDINATES in Fc
m2 := <Xi[2],Yi[2],Zi[2]>:
m3 := <Xi[3],Yi[3],Zi[3]>:
m4 := <Xi[4],Yi[4],Zi[4]>:
#=================================================================================

#=================================================================================
# DISTANCE CONSTRAINTS:
# The distance between each two points expressed in the (s,Z)-space
c1 := Transpose(m2-m1).(m2-m1) - d12:
c2 := Transpose(m3-m1).(m3-m1) - d13:
c3 := Transpose(m4-m1).(m4-m1) - d14:
c4 := Transpose(m3-m2).(m3-m2) - d23:
c5 := Transpose(m4-m2).(m4-m2) - d24:
c6 := Transpose(m4-m3).(m4-m3) - d34:
const := map(expand, [c1,c2,c3,c4,c5,c6]):
const_ := [c1,c2,c3,c4,c5,c6]:
#=================================================================================

#=================================================================================
# SYSTEM OF EQUATIONS in (s,Z)-space
err := Vector(2*k):   # ERROR VECTOR
for i from 1 to k do
  err[2*i-1] := xi[i] - xi_star[i]:
  err[2*i] := yi[i] - yi_star[i]:
od:
err := map(normal, err):

Lp := Matrix(2*k, 6):   # INTERACTION MATRIX
for i from 1 to k do
  Lp_i := <-1/Zi[i], 0, xi[i]/Zi[i], xi[i]*yi[i], -(1 + xi[i]^2), yi[i]>:
  Lp_i_p1 := <0, -1/Zi[i], yi[i]/Zi[i], (1 + yi[i]^2), -xi[i]*yi[i], -xi[i]>:
  Lp(2*i-1, 1..) := Lp_i:
  Lp(2*i, 1..) := Lp_i_p1:
od:

F := map(expand, map(numer, map(normal, convert(Transpose(err).Lp, list)))):    # GRADIENT of V=err^T.err
#=================================================================================

#=================================================================================
# COPLANARITY CONDITION

# J1 is the determinant of matrix
# | X1 X2 X3 X4 |
# | Y1 Y2 Y3 Y4 |
# | Z1 Z2 Z3 Z4 |
# |  1  1  1  1 |
# J1=0 iff the points are coplanar: it can then be used to accelerate the computations

J1 := Z1*Z2*Z3*x1*y2-Z1*Z2*Z3*x1*y3-Z1*Z2*Z3*x2*y1+Z1*Z2*Z3*x2*y3+Z1*Z2*Z3*x3*y1-Z1*Z2*Z3*x3*y2-
Z1*Z2*Z4*x1*y2+Z1*Z2*Z4*x1*y4+Z1*Z2*Z4*x2*y1-Z1*Z2*Z4*x2*y4-Z1*Z2*Z4*x4*y1+Z1*Z2*Z4*x4*y2+
Z1*Z3*Z4*x1*y3-Z1*Z3*Z4*x1*y4-Z1*Z3*Z4*x3*y1+Z1*Z3*Z4*x3*y4+Z1*Z3*Z4*x4*y1-Z1*Z3*Z4*x4*y3-
Z2*Z3*Z4*x2*y3+Z2*Z3*Z4*x2*y4+Z2*Z3*Z4*x3*y2-Z2*Z3*Z4*x3*y4-Z2*Z3*Z4*x4*y2+Z2*Z3*Z4*x4*y3:

(*
# --------------------------------------------------------------------------------
# TEST that J1 divides the determinant of Jac(sys)
# The factorization step is quite long

sys := [op(F), op(const)]:
vars := [op(indets(sys))]:

jac_sys := Matrix(12,12):
for i from 1 to numelems(sys) do
  jac_sys[i,1] := diff(sys[i], x1):
  jac_sys[i,2] := diff(sys[i], x2):
  jac_sys[i,3] := diff(sys[i], x3):
  jac_sys[i,4] := diff(sys[i], x4):
  jac_sys[i,5] := diff(sys[i], y1):
  jac_sys[i,6] := diff(sys[i], y2):
  jac_sys[i,7] := diff(sys[i], y3):
  jac_sys[i,8] := diff(sys[i], y4):
  jac_sys[i,9] := diff(sys[i], Z1):
  jac_sys[i,10] := diff(sys[i], Z2):
  jac_sys[i,11] := diff(sys[i], Z3):
  jac_sys[i,12] := diff(sys[i], Z4):
od:

jac_sys_det := Determinant(jac_sys):

det_factors := map(l->l[1],factors(jac_sys_det)[2]):
min_factor := op(map(l-> if degree(l)=5 then l: fi, det_factors )):

# Compare min_factor (factor of smallest degree) with J1
# --------------------------------------------------------------------------------
*)

#=================================================================================
# COMPUTATION OF THE CRITICAL POINTS

# generate a RANDOM LINEAR FORM - introduces genericity which helps Grobner Bases computations
randomize();
rand_pol := randpoly([t, l, op(indets(F))], degree=1, dense, homogeneous):

# SYSTEM OF EQUATIONS
#sys := [ rand_pol, 1-l*Z1*Z2*Z3*Z4, op(F), op(const)]: # without COPLANARITY condition
sys := [ rand_pol, 1-l*Z1*Z2*Z3*Z4, J1, op(F), op(const)]: # WITH COPLANARITY condition

vars := [t, l, op(indets(F))]:    #VARIABLES

# SOLVING the system with msolve:
ToMSolve(sys, 0, vars, fname1);  # generate msolve input file

# EITHER run the next line in the terminal:
#/home/polsys2/garciaf/msolve -v 2 -t 12 -f /tmp/jorge.rectangle0_sZ.sys -o /tmp/jorge.rectangle0_sZ.sols

# OR exectute in background of current maple session with:
str := cat("/home/polsys2/garciaf/msolve.experimental -v 2 -t 12 -f ", fname1, " -o ", fname2, " > ", fname_log, " 2>&1 &");
system(str):

# The file fname_log will contain information such as:
# degree of the ideal, number of complex & real solutions, computing time
# The file fname2 will contain a zero-dimensional rational parametrization of the solution set - encoded in msolve-format
# go to https://msolve.lip6.fr/ to learn more about the rational parametrization of the zero-dim variety
#=================================================================================

#=================================================================================
# PROCESSING SOLUTIONS
#=================================================================================

#=================================================================================
# Retrieve SOLUTIONS from msolve OUTPUT

fname_output := fname2:
read(fname_output):
sols:=map(s->[seq(vars[i]=s[i],i=1..nops(vars))], %[2][2]):   # vector of all REAL SOLUTIONS of the system
vars := [t, l, op(indets(F))]:

(*
# --------------------------------------------------------------------------------
# VERIFICATION step: evaluate polynomials at the computed solutions - the result should be very small, else something failed
evalf(subs( map(l-> lhs(l) = rhs(l)[1], sols[1]) , sys));   # evaluate sys at the FIRST SOLUTION ONLY

# absolute MAX of all the polynomials (excluding the random linear form) EVALUATED at all solutions -
max_eval := 0:
for i from 1 to numelems(sols) do
  max_eval := max([max_eval, max(map(abs, evalf(subs( map(l-> lhs(l) = rhs(l)[1], sols[i]) , sys[2..-1]))))]);
od:
max_eval;   # if max_eval is large, smth went WRONG
# --------------------------------------------------------------------------------
*)

#=================================================================================
# REMOVE solutions BEHIND IMAGE PLANE
# i.e. Keep only the points that are visible (Zi>0 for all i)

Z_coords := []:
Z_visible := []:
inds := []:
for i from 1 to numelems(sols) do
  Z_coords := [op(Z_coords), map( l->rhs(l)[1] , sols[i][3..6])]:
  flagg := true;
  for ii from 1 to 4 do   # for each solution check that all Zi > 0
      if Z_coords[i][ii]<0 then flagg:=false: fi:
  od:
  if flagg then inds := [op(inds), i]: fi:
od:
Z_visible := [seq(Z_coords[i], i = inds)]:

sols_vis := []:   # keep only VISIBLE (real) SOLUTIONS
for i in inds do
sols_vis := [op(sols_vis), map(l->lhs(l)=rhs(l)[1], sols[i][1..])]:
od:
sols_vis_ := evalf(sols_vis):

s_sols := map( l-> map(rhs, l[7..14]) , sols_vis):    # s-coordinates of the solutions
Z_sols := map( l-> map(rhs, l[3..6]) , sols_vis):     # Z-coordinates
#=================================================================================

#=================================================================================
# RETRIEVING CAMERA POSE (t, R): position vector+rotation matrix
# need to differentiate between COPLANAR and NON-COPLANAR objects
interface(displayprecision = 5);

# --------------------------------------------------------------------------------
# PLANAR objects
Plin := < OMi[1][1],OMi[1][2], 0,0, 0,0, 1,0,0;
          0,0, OMi[1][1],OMi[1][2], 0,0, 0,1,0;
          0,0, 0,0, OMi[1][1],OMi[1][2], 0,0,1;
          OMi[2][1],OMi[2][2], 0,0, 0,0, 1,0,0;
          0,0, OMi[2][1],OMi[2][2], 0,0, 0,1,0;
          0,0, 0,0, OMi[2][1],OMi[2][2], 0,0,1;
          OMi[3][1],OMi[3][2], 0,0, 0,0, 1,0,0;
          0,0, OMi[3][1],OMi[3][2], 0,0, 0,1,0;
          0,0, 0,0, OMi[3][1],OMi[3][2], 0,0,1;
          OMi[4][1],OMi[4][2], 0,0, 0,0, 1,0,0;
          0,0, OMi[4][1],OMi[4][2], 0,0, 0,1,0;
          0,0, 0,0, OMi[4][1],OMi[4][2], 0,0,1>:

i:= 1:  # SELECT SOLUTION NUMBER
# For each solution, we retrieve the camera pose and classify it in minima, maxima or local minima

Xlin := evalf(< s_sols[i][1]*Z_sols[i][1], s_sols[i][5]*Z_sols[i][1], Z_sols[i][1],
        s_sols[i][2]*Z_sols[i][2], s_sols[i][6]*Z_sols[i][2], Z_sols[i][2],
        s_sols[i][3]*Z_sols[i][3], s_sols[i][7]*Z_sols[i][3], Z_sols[i][3],
        s_sols[i][4]*Z_sols[i][4], s_sols[i][8]*Z_sols[i][4], Z_sols[i][4] >):

# solve REDUCED 3Nx9 system
That_sol := LinearSolve(Plin, Xlin):
t_sol := evalf(That_sol[7..9]):
R_red := evalf(Transpose(<That_sol[1..2] | That_sol[3..4] | That_sol[5..6]>)):  # first two columns of Rot matrix

# To each solution correspond 2 points in E(3)
# one is a true camera pose and one is its reflection about the plane of the object

sysR := convert(<Transpose(R_red).<r13_,r23_,r33_>; Transpose(<r13_,r23_,r33_>).<r13_,r23_,r33_> -1> , list);
r_candidates := solve(sysR, [r13_,r23_,r33_]):  # solve quadratic eqs to retrieve the last column of Rot matrix

# r_candidates contains 2 solutions - need to select that which makes detR = 1

r_3row := map(rhs, r_candidates[2]);  # CHANGE for r_candidates[2] IF detR=-1
r13 := r_3row[1]: r23 := r_3row[2]: r33 := r_3row[3]:
R_sol := evalf(Transpose(<That_sol[1..2], r13 | That_sol[3..4], r23 | That_sol[5..6], r33>)):
detR := LinearAlgebra:-Determinant(R_sol);
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# NON-PLANAR OBJECTS
Plin := < OMi[1][1],OMi[1][2],OMi[1][3], 0,0,0, 0,0,0, 1,0,0;
          0,0,0, OMi[1][1],OMi[1][2],OMi[1][3], 0,0,0, 0,1,0;
          0,0,0, 0,0,0, OMi[1][1],OMi[1][2],OMi[1][3], 0,0,1;
          OMi[2][1],OMi[2][2],OMi[2][3], 0,0,0, 0,0,0, 1,0,0;
          0,0,0, OMi[2][1],OMi[2][2],OMi[2][3], 0,0,0, 0,1,0;
          0,0,0, 0,0,0, OMi[2][1],OMi[2][2],OMi[2][3], 0,0,1;
          OMi[3][1],OMi[3][2],OMi[3][3], 0,0,0, 0,0,0, 1,0,0;
          0,0,0, OMi[3][1],OMi[3][2],OMi[3][3], 0,0,0, 0,1,0;
          0,0,0, 0,0,0, OMi[3][1],OMi[3][2],OMi[3][3], 0,0,1;
          OMi[4][1],OMi[4][2],OMi[4][3], 0,0,0, 0,0,0, 1,0,0;
          0,0,0, OMi[4][1],OMi[4][2],OMi[4][3], 0,0,0, 0,1,0;
          0,0,0, 0,0,0, OMi[4][1],OMi[4][2],OMi[4][3], 0,0,1>:

i:= 1:  # SELECT SOLUTION NUMBER

Xlin := evalf(< s_sols[i][1]*Z_sols[i][1], s_sols[i][5]*Z_sols[i][1], Z_sols[i][1],
          s_sols[i][2]*Z_sols[i][2], s_sols[i][6]*Z_sols[i][2], Z_sols[i][2],
          s_sols[i][3]*Z_sols[i][3], s_sols[i][7]*Z_sols[i][3], Z_sols[i][3],
          s_sols[i][4]*Z_sols[i][4], s_sols[i][8]*Z_sols[i][4], Z_sols[i][4] >):

# solve FULL 3Nx12 system
Mlin_sol := LinearAlgebra:-MatrixInverse(Plin).Xlin:

R_sol := evalf(Transpose(<Mlin_sol[1..3] | Mlin_sol[4..6] | Mlin_sol[7..9]>)):
t_sol := evalf(Mlin_sol[10..12]):

# Some solutions will correspond to a true camera pose in SE(3)
# and some will correspond to a reflected pose in E^-(3) (indirect Euclidean isometries)
# verify by checking that detR=1

detR := LinearAlgebra:-Determinant(R_sol);  # if detR=-1, DISCARD solution
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# PRINT solutions on screen
evalf(<sd[1..4];sd[5..8]>);
s_table := evalf[5](<s_sols[i][1..4];s_sols[i][5..8]>);
Z_table := evalf[5](Z_sols[i]);
t_table := map(Round, t_sol, 5);
R_table := map(Round, R_sol, 5);
# --------------------------------------------------------------------------------

#=================================================================================
# CLASSIFYING SOLUTIONS IN MINIMA, MAXIMA OR SADDLE POINTS
#=================================================================================

#=================================================================================
# COMPUTE tu-vector corresponding to Rot matrix

# Inverse Rodrigues Formula
eigen_R := Eigenvectors(R_sol);
u_axis := Re(Column(eigen_R[2], 1)); # Select eigenvector corresponding to eigenvalue 1
arccos( (Trace(R_sol)-1)/2 );
theta := Re(%); # may need a sign change
#theta := Re(eigen_R[1][1]);

R_sol;  # VERIFY that both are equal, otherwise change sign of theta
Rot_test := Student[LinearAlgebra]:-RotationMatrix(theta, u_axis);
#=================================================================================

#=================================================================================
# Compute constrained HESSIAN of V
# Need to express V as a function of t=(x,y,z) and q=(r,u,v,w)

Rot_quat := < r^2+u^2-v^2-w^2,  2*(u*v-r*w),     2*(r*v+u*w) ; # ROTATION matrix expressed in QUATERNION components
            2*(r*w+u*v),     r^2-u^2+v^2-w^2,  2*(v*w-r*u) ;
            2*(u*w-r*v),     2*(r*u+v*w),     r^2-u^2-v^2+w^2> :
const_quat := r^2+u^2+v^2+w^2 - 1:    # unit quaternion constraint
tc := <x,y,z>:    # vector CO in camera frame
vars := [r, u, v, w, x, y, z]:

CMi := map(l->(Rot_quat.l + tc), OMi ):
Xi := map(l->l[1], CMi): Yi := map(l->l[2], CMi): Zi := map(l->l[3], CMi):
xi := Xi/~Zi: yi := Yi/~Zi:

err := [ op(xi-~sd[1..4]), op(yi-~sd[5..8]) ]:  # error vector
V := normal(add(err *~ err)):

dV_dr := diff(V, r): dV_du := diff(V, u): dV_dv := diff(V, v): dV_dw := diff(V, w):
dV_dx := diff(V, x): dV_dy := diff(V, y): dV_dz := diff(V, z):

grad_V := normal([dV_dr, dV_du, dV_dv, dV_dw, dV_dx, dV_dy, dV_dz]):  # GRADIENT with respect to (t,q)

(*
# HESSIAN with respect to (t,q) - unconstrained (7x7)-Hessian
# Computing the Hessian symbolically is expensive - it should be done only once, not for every selected solution
dV2_drdr := diff(dV_dr, r): dV2_drdu := diff(dV_dr, u): dV2_drdv := diff(dV_dr, v): dV2_drdw := diff(dV_dr, w): dV2_drdx := diff(dV_dr, x): dV2_drdy := diff(dV_dr, y): dV2_drdz := diff(dV_dr, z):
dV2_dudu := diff(dV_du, u): dV2_dudv := diff(dV_du, v): dV2_dudw := diff(dV_du, w): dV2_dudx := diff(dV_du, x): dV2_dudy := diff(dV_du, y): dV2_dudz := diff(dV_du, z):
dV2_dvdv := diff(dV_dv, v): dV2_dvdw := diff(dV_dv, w): dV2_dvdx := diff(dV_dv, x): dV2_dvdy := diff(dV_dv, y): dV2_dvdz := diff(dV_dv, z):
dV2_dwdw := diff(dV_dw, w): dV2_dwdx := diff(dV_dw, x): dV2_dwdy := diff(dV_dw, y): dV2_dwdz := diff(dV_dw, z):
dV2_dxdx := diff(dV_dx, x): dV2_dxdy := diff(dV_dx, y): dV2_dxdz := diff(dV_dx, z):
dV2_dydy := diff(dV_dy, y): dV2_dydz := diff(dV_dy, z):
dV2_dzdz := diff(dV_dz, z):

hess_V := map(normal,<
dV2_drdr, dV2_drdu, dV2_drdv, dV2_drdw, dV2_drdx, dV2_drdy, dV2_drdz ;
dV2_drdu, dV2_dudu, dV2_dudv, dV2_dudw, dV2_dudx, dV2_dudy, dV2_dudz ;
dV2_drdv, dV2_dudv, dV2_dvdv, dV2_dvdw, dV2_dvdx, dV2_dvdy, dV2_dvdz ;
dV2_drdw, dV2_dudw, dV2_dvdw, dV2_dwdw, dV2_dwdx, dV2_dwdy, dV2_dwdz ;
dV2_drdx, dV2_dudx, dV2_dvdx, dV2_dwdx, dV2_dxdx, dV2_dxdy, dV2_dxdz ;
dV2_drdy, dV2_dudy, dV2_dvdy, dV2_dwdy, dV2_dxdy, dV2_dydy, dV2_dydz ;
dV2_drdz, dV2_dudz, dV2_dvdz, dV2_dwdz, dV2_dxdz, dV2_dydz, dV2_dzdz >):
#  hess_Vn := map(numer, hess_V):
*)

# MATRIX FOR CHANGE OF BASIS FOR THE HESSIAN computing the constrained (6x6)-Hessian
c := r^2+u^2+v^2+w^2 - 1: # unit quaternion constraint
dc_dr := diff(c, r): dc_du := diff(c, u): dc_dv := diff(c, v): dc_dw := diff(c, w): dc_dx := diff(c, x): dc_dy := diff(c, y): dc_dz := diff(c, z):
grad_c := Matrix(LinearAlgebra:-Transpose(< dc_dr, dc_du, dc_dv, dc_dw, dc_dx, dc_dy, dc_dz >)):
nullspaceC := LinearAlgebra:-NullSpace(grad_c):

Zc := Matrix(7,6):
for j from 1 to numelems(nullspaceC) do
  Zc[1..,j] := nullspaceC[j]:
od:
#=================================================================================

#===============================================================================
# EVALUATE ERROR, and the EIGENVALUES of HESS(V) at the SELECTED SOLUTION
# It must be repeated for every solution to classify the local minima and saddle points

err_crit := eval( err, [r = cos(theta/2), u = u_axis[1]*sin(theta/2),
v = u_axis[2]*sin(theta/2), w = u_axis[3]*sin(theta/2),
x = t_sol[1], y = t_sol[2], z = t_sol[3]] );

hess_Vcrit  := eval( hess_V, [r = cos(theta/2), u = u_axis[1]*sin(theta/2),
v = u_axis[2]*sin(theta/2), w = u_axis[3]*sin(theta/2),
x = t_sol[1], y = t_sol[2], z = t_sol[3]] ):

Zccrit  := eval( Zc, [r = cos(theta/2), u = u_axis[1]*sin(theta/2),
v = u_axis[2]*sin(theta/2), w = u_axis[3]*sin(theta/2),
x = t_sol[1], y = t_sol[2], z = t_sol[3]] ):

H_V_crit := Transpose(Zccrit).hess_Vcrit.Zccrit:
Eigenvalues(H_V_crit);  # EIGENVALUES of (constrained) Hessian
# If the eigenvalues are all positive, the point is a local minimum
# if at least one eigenvalue is negative, it is a saddle point (there are no local maxima)

Round:= (x,n)-> parse(sprintf(sprintf("%%.%df",n),x));
eigenvalues := map(Round, convert(Re(Eigenvalues(H_V)), list) , 4);
#===============================================================================

#===============================================================================
# THE END
#===============================================================================
