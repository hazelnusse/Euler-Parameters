autoz on
autorhs on
newtonian N
body A
variable e0', e{3}'
motionvariable' u0', u{2}'
constants ma, Ixx, Iyy, Izz, Ixy, Iyz, Ixz
constants Tax, Tay, Taz

% Set mass and inertia of rigidy body A
mass a=ma
inertia a, Ixx, Iyy, Izz, Ixy, Iyz, Ixz

% Orient A relative to N
dircos(N, A, euler, e0, e1, e2, e3)

% Declare angular velocity of A relative to N
w_a_n> = u0*a1> + u1*a2> + u2*a3>

% Ignore translational motion
v_ao_n> = 0>
a_ao_n> = 0>

% Form kinematic equations associated with orientation
kindiffs(N, A, euler, e0, e1, e2, e3)

% Form angular acceleration of A in N
alf_a_n> = dt(w_a_n>, n)

% Declare applied torques acting on A
torque_a> = Tax*a1> + Tay*a2> + Taz*a3>

% Form Kane's dynamical equations
zero = fr() + frstar()

% Solve them for the time derivatives of the generalized speeds
solve(zero, [u0', u1', u2'])

% Form that A and B matrices by computing partial derivatives
states = [e0, e1, e2, e3, u0, u1, u2]
f = [rhs(e0'); rhs(e1'); rhs(e2'); rhs(e3'); rhs(u0'); rhs(u1'); rhs(u2')]
A = d(f, states)
B = [0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; d(f[5], Tax), d(f[5], Tay), d(f[5], Taz); d(f[6], Tax), d(f[6], Tay), d(f[6], Taz); d(f[7], Tax), d(f[7], Tay), d(f[7], Taz)]

% 4x4 transformation matrix used by OpenGL, column-major format
m = [N_A[1,1], N_A[2,1], N_A[3,1], 0, N_A[1,2], N_A[2,2], N_A[3,2], 0, N_A[1,3], N_A[2,3], N_A[3,3], 0, 0, 0, 0, 1]

unitsystem kg, m, s
input ma=1, Ixx=1, Iyy=1, Izz=1, Ixy=1, Iyz=1, Ixz=1
input Tax=0, Tay=0, Taz=0
input e0=0, e1=0, e2=0, e3=1, u0=0, u1=0, u2=0

output e0, e1, e2, e3, u0, u1, u2, m, A, B
code dynamics() rigidbody.c
