function [u_opt, xk_opt, xkp_opt, sol] = decoupledQP(k,zK,zKP,lambdaK,lambdaKP,mpc_problem, optimization_options, aladin_options)
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

nx = size(A,1); 
nu = size(B,2); 
%% Defining and solving optimization problem
x_k = sdpvar(nx,1);
x_kp = sdpvar(nx,1);
u_k = sdpvar(nu,1);

if (k == 1)
    objective = x_k'*Q*x_k + u_k'*R*u_k - lambdaKP'*(x_kp-zKP) + rho/2*(x_kp-zKP)'*omega*(x_kp-zKP);  
    constraints = [x_kp == A*x_k + B*u_k, umin <= u_k <= umax, xmin <= x_k <= xmax, xmin <= x_kp <= xmax; x_k == x0];
elseif ((k < N) & (k > 1)) 
    objective = x_k'*Q*x_k + u_k'*R*u_k + lambdaK'*(zK-x_k) - lambdaKP'*(x_kp-zKP)+...
        rho/2*(x_k - zK)'*omega*(x_k - zK)+rho/2*(x_kp-zKP)'*omega*(x_kp-zKP);  
    constraints = [x_kp == A*x_k + B*u_k, umin <= u_k <= umax, xmin <= x_k <= xmax, xmin <= x_kp <= xmax ];
elseif (k == N)
    objective = x_k'*Q*x_k + u_k'*R*u_k + lambdaK'*(zK-x_k) + rho/2*(x_k - zK)'*omega*(x_k - zK);  
    constraints = [x_kp == A*x_k + B*u_k, umin <= u_k <= umax, xmin <= x_k <= xmax, xmin <= x_kp <= xmax ];
else
    error(sprintf('DecoupledQP: Unexpected step k = %d', k));
end

sol = optimize(constraints,objective,optim_setup);

u_opt = value(u_k);
xk_opt =value(x_k);
xkp_opt =value(x_kp);
end