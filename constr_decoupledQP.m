function [mpsol1, mpsol2, mpsol3, DIAGNOSTIC1, DIAGNOSTIC2, DIAGNOSTIC3] = constr_decoupledQP(mpc_problem, optimization_options, aladin_options)

%% Loading system specifications

A = mpc_problem.A;
B = mpc_problem.B;
umin = mpc_problem.umin;
umax = mpc_problem.umax;
xmin = mpc_problem.xmin;
xmax = mpc_problem.xmax;
x0 = mpc_problem.x0;

%% Loading optimization settings
Q = mpc_problem.Q;
R = mpc_problem.R;
N = mpc_problem.N;
rho = aladin_options.rho;
omega = aladin_options.omega;

yalmip('clear');
my_solver = optimization_options.yalmip.solver;
my_verbose = optimization_options.yalmip.verbose;
optim_setup = sdpsettings('solver',my_solver,'verbose',my_verbose);

rho = aladin_options.rho;
maxiter = aladin_options.maxiter;
init_var = aladin_options.init_var;
tolerance = aladin_options.tolerance;

nx = size(A,1); 
nu = size(B,2); 

%% Defining and solving optimization problem

x_k = sdpvar(nx,1);
x_kp = sdpvar(nx,1);
u_k = sdpvar(nu,1);

zK = sdpvar(nx,1);
zKP = sdpvar(nx,1);
lambdaK = sdpvar(nx,1);
lambdaKP = sdpvar(nx,1);

%% Initial step decoupled QP

objective = x_k'*Q*x_k + u_k'*R*u_k - lambdaKP'*(x_kp-zKP) + rho/2*(x_kp-zKP)'*omega*(x_kp-zKP);  
constraints = [x_kp == A*x_k + B*u_k, umin <= u_k <= umax, xmin <= x_k <= xmax, xmin <= x_kp <= xmax, x_k == zK];
constraints = [constraints, xmin <= zKP <= xmax, -10 <= lambdaKP <= 1e6];

[mpsol1, DIAGNOSTIC1,aux1,valueFunction1,ZPWF1] = solvemp(constraints,objective,[], [zK; zKP; lambdaKP], [u_k; x_k; x_kp]);
[ mpsol1 ] = get_explicit_partition( mpsol1 );

%% Middle steps decoupled QP

objective = x_k'*Q*x_k + u_k'*R*u_k + lambdaK'*(zK-x_k) - lambdaKP'*(x_kp-zKP)+...
    rho/2*(x_k - zK)'*omega*(x_k - zK)+rho/2*(x_kp-zKP)'*omega*(x_kp-zKP);  
constraints = [x_kp == A*x_k + B*u_k, umin <= u_k <= umax, xmin <= x_k <= xmax, xmin <= x_kp <= xmax ];
constraints = [constraints, xmin <= zK <= xmax, xmin <= zKP <= xmax, -10 <= lambdaK <= 1e6, -10 <= lambdaK <= 1e6];
[mpsol2, DIAGNOSTIC2,aux2,valueFunction2,ZPWF2] = solvemp(constraints,objective,[], [zK; zKP; lambdaK; lambdaKP], [u_k; x_k; x_kp]);
[ mpsol2 ] = get_explicit_partition( mpsol2 );

%% Final step decoupled QP

objective = x_k'*Q*x_k + u_k'*R*u_k + lambdaK'*(zK-x_k) + rho/2*(x_k - zK)'*omega*(x_k - zK);  
constraints = [x_kp == A*x_k + B*u_k, umin <= u_k <= umax, xmin <= x_k <= xmax, xmin <= x_kp <= xmax ];
constraints = [constraints, xmin <= zK <= xmax, -10 <= lambdaK <= 1e6];
[mpsol3, DIAGNOSTIC3,aux3,valueFunction3,ZPWF3] = solvemp(constraints,objective,[], [zK; lambdaK], [u_k; x_k; x_kp]);
[ mpsol3 ] = get_explicit_partition( mpsol3 );

end