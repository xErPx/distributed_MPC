function solution = ModelPredictiveControl(mpc_problem, optimization_options, solution)
% function [u_opt, x_opt, sol_MPC] = ModelPredictiveControl(mpc_problem, optimization_options)

linewidth = 1;
markersize = 6;
mpccolor = '#808080';
% fontsize = 30;

%% Loading system specifications
A = mpc_problem.A;
B = mpc_problem.B;
umin = mpc_problem.umin;
umax = mpc_problem.umax;
xmin = mpc_problem.xmin;
xmax = mpc_problem.xmax;
x0 = mpc_problem.x0;
xs = [mpc_problem.x1s; mpc_problem.x2s];
us = mpc_problem.us;

%% Loading optimization settings
Q = mpc_problem.Q;
R = mpc_problem.R;
N = mpc_problem.N;
aladin_iter = solution.sumiter;

yalmip('clear');
mySolver = optimization_options.yalmip.solver;
myVerbose = optimization_options.yalmip.verbose;
graph = optimization_options.graph;
optim_options = sdpsettings('solver',mySolver,'verbose',myVerbose);

nx = size(A,1);
nu = size(B,2);

%% Defining and solving optimization problem
u = cell(N,1);          
x = cell(N+1,1);     

for k = 1:N+1
    x{k} = sdpvar(nx,1);    
    if k <= N
        u{k} = sdpvar(nu,1);
    end
end

constraints = [x{1} == x0];
objective = 0;

for k = 1:N
    objective = objective + (x{k}'*Q*x{k}) + u{k}'*R*u{k};  
    constraints = [constraints, x{k+1} == A*x{k}+B*u{k}];
    constraints = [constraints, umin<=u{k}<=umax, xmin <= x{k} <= xmax]; 
end

sol_MPC  = optimize(constraints, objective, optim_options);
for k =1:N+1
    x_opt{k} = value(x{k});
    if k <= N
        u_opt{k} = value(u{k});
    end
end

solution.x_mpc = x_opt;
solution.u_mpc = u_opt;
solution.sol_mpc = sol_MPC;

if graph == 1
    %% Graphs

    path = [cd, '\Figures\plot'];
    
    % state variables graph
    for k = 1:N+1
        x_opt_vector = cell2mat(x_opt);
    end
    for j = 1:aladin_iter
        figure(j)
        for k = 1:1:nx
            subplot(nx, 2, k)
            hold on;
            box on; 
            grid on;
            stairs([0:N],x_opt_vector(k,:)+xs(k),'Color', mpccolor, 'LineStyle', '--','LineWidth',linewidth,'DisplayName','MPC solution');
            axis padded;
            xlabel('prediction horizon');ylabel(['$h_' num2str(k) '$ [m]']);
        end
    end
    
    % Input variable graph
    for k = 1:N
        u_opt_vector = cell2mat(u_opt);
    end
    
    for j = 1:aladin_iter
        figure(j)
        subplot(nx,2,nx+1)
        hold on; box on; grid on;
        a = stairs([0:N-1],[u_opt_vector(1,:)+us],'Color', mpccolor, 'LineStyle', '--','LineWidth',linewidth,'DisplayName','MPC solution');
        uistack(a,"bottom");
        axis padded;
        xlabel('prediction horizon');ylabel('$q_0$ [m$^3$/s]');
    end
end