autoz off
newtonian N
body A
variable e0', e{3}', x', y', z'
motionvariable' u0', u{5}'
constants ma, Ixx, Iyy, Izz, Ixz
constants g, Fax, Fay, Faz, Fnx, Fny, Fnz, Tax, Tay, Taz, Tnx, Tny, Tnz

% Set mass and inertia of rigidy body A
mass a=ma
inertia a, Ixx, Iyy, Izz, 0, 0, Ixz

% Orient A relative to N
dircos(N, A, euler, e0, e1, e2, e3)

% Declare angular velocity of A relative to N
w_a_n> = u0*a1> + u1*a2> + u2*a3>

% Declare velocity of mass center of A as viewed by an observer fixed in N
v_ao_n> = u3*a1> + u4*a2> + u5*a3>

% Form kinematic equations associated with orientation
kindiffs(N, A, euler, e0, e1, e2, e3)

% Form kinematic equations associated with translation
x' = dot(v_ao_n>, n1>)
y' = dot(v_ao_n>, n2>)
z' = dot(v_ao_n>, n3>)

% Form angular acceleration of A in N
alf_a_n> = dt(w_a_n>, n)

% Form acceleration of mass center of A as viewed by an observer fixed in N
a_ao_n> = dt(v_ao_n>, n)

% Declare applied torques acting on A
torque_a> = Tax*a1> + Tay*a2> + Taz*a3> + Tnx*n1> + Tny*n2> + Tnz*n3>

% Declare applied forces acting at mass center of A
force_ao> = Fax*a1> + Fay*a2> + Faz*a3> + Fnx*n1> + Fny*n2> + Fnz*n3>

% Form Kane's dynamical equations
zero = fr() + frstar()

% Solve them for the time derivatives of the generalized speeds
solve(zero, [u0', u1', u2', u3', u4', u5'])

% 4x4 transformation matrix used by OpenGL, column-major format
m = [N_A[1,1], N_A[2,1], N_A[3,1], 0, N_A[1,2], N_A[2,2], N_A[3,2], 0, N_A[1,3], N_A[2,3], N_A[3,3], 0, x, y, z, 1]

ke = (ma*dot(v_ao_n>, v_ao_n>) + dot(w_a_n>, dot(I_A_AO>>, w_a_n>)))/2.0
pe = ma*g*z
te = ke + pe

unitsystem kg, m, s
input ma=1,Ixx=1,Iyy=1,Izz=1,Ixz=1
input Fax=0,Fay=0,Faz=0,Fnx=0,Fny=0,Fnz=0,Tax=0,Tay=0,Taz=0,Tnx=0,Tny=0,Tnz=0
input e0=0, e1=0, e2=0, e3=1, x=0, y=0, z=0, u0=0, u1=0, u2=0, u3=0, u4=0, u5=0

output e0, e1, e2, e3, x, y, z, u0, u1, u2, u3, u4, u5, m, ke, pe, te
code dynamics() rigidbody.c
