function [dz_opt, du_opt, lambda_opt, sol] = coupledQP(xkp_opt, u_opt, mpc_problem,optimization_options,aladin_options)
%% Loading system specifications
A = mpc_problem.A;
B = mpc_problem.B;
xp = xkp_opt;

%% Loading optimization settings
Q = mpc_problem.Q;
R = mpc_problem.R;
N = mpc_problem.N;

yalmip('clear');
my_solver = optimization_options.yalmip.solver;
my_verbose = optimization_options.yalmip.verbose;
optim_setup = sdpsettings('solver',my_solver,'verbose',my_verbose);

nx = size(A,1); 
nu = size(B,2); 

%% Defining and solving optimization problem
dz = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
du = sdpvar(repmat(nu,1,N),repmat(1,1,N));

% Zeros = zeros(length(nx),1);
objective = 0;
constraints = [dz{1} == 0];
for i=1:N
    objective = objective + dz{i}'*Q*dz{i} + du{i}'*R*du{i};
    constraints = [constraints, xp{i+1} + dz{i+1} == A*(xp{i} + dz{i}) + B*(u_opt{i} + du{i})];
end

sol = optimize(constraints,objective,optim_setup);

for k =1:length(du)
    dz_opt{k} = value(dz{k});
    du_opt{k} = value(du{k});
    lambda_opt{k} = -dual(constraints(k));
end
dz_opt{k+1} = value(dz{k+1});
lambda_opt{k+1} = -dual(constraints(k+1));
end