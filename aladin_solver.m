function [solution] = aladin_solver(mpc_problem, optimization_options, aladin_options)

%% Visualisation settings
linewidth = 1;
markersize = 6;
impcolor = 'blue';

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
tolerance = aladin_options.tolerance;
maxiter = aladin_options.maxiter;
rho = aladin_options.rho;
init_var = aladin_options.init_var;
graph = optimization_options.graph;

nx = size(A,1); 
nu = size(B,2); 

%% Inicialization of algorithm parameters and variables using LQR or zeros

[K,P] = dlqr(A,B,Q,R);
if(isequal (init_var,'dlqr'))
    z{1}{1} = x0;
    for k=1:N
        u{k} = -K*z{1}{k};
        z{1}{k+1}=A*z{1}{k}+B*u{k};
    end
elseif(isequal (init_var,'zeros'))
    z{1}{1} = x0;
    for i = 1:N+1
        z{1}{i}=zeros(nx,1);
    end
else
    error(sprintf('Unexpected input: %s!', init_var))
end

omega = P; %P -> Riccati equation solution
aladin_options.omega = omega;

for i=1:N+1
    lambda{1}{i} = [0; 0];
end

j = 1;
dz{1}{1} = ones(nx,1)*Inf;

%% Construction of explicit MPC

if aladin_options.sol_method == 1

    [mpsol1, mpsol2, mpsol3, diagnostic1, diagnostic2, diagnostic3] = constr_decoupledQP(mpc_problem, optimization_options, aladin_options);
    solution.mpsol1 = mpsol1;
    solution.mpsol2 = mpsol2;
    solution.mpsol3 = mpsol3;
    % 
    solution.diagnostic{1} = diagnostic1;
    solution.diagnostic{2} = diagnostic2;
    solution.diagnostic{3} = diagnostic3;
    
end

%% ALADIN algorithm iterations

while (j <= maxiter) & (rho*norm(cell2mat(dz{j}),1) > tolerance)
    xkp_opt{j}{1} = x0;
    
    for i = 1:N
        if aladin_options.sol_method == 0
            [u_opt{j}{i}, xk_opt{j}{i}, xkp_opt{j}{i+1}, sol_decoupled{j}{i}] = decoupledQP(i,z{j}{i},z{j}{i+1},lambda{j}{i},lambda{j}{i+1},mpc_problem,optimization_options,aladin_options);
        elseif aladin_options.sol_method == 1
            [u_opt{j}{i}, xk_opt{j}{i}, xkp_opt{j}{i+1}, sol_decoupled{j}{i}] = decoupled_MPQP(i,z{j}{i},z{j}{i+1},lambda{j}{i},lambda{j}{i+1},solution, mpc_problem);
        elseif aladin_options.sol_method == 2
            [u_opt{j}{i}, xk_opt{j}{i}, xkp_opt{j}{i+1}, sol_decoupled{j}{i}] = decoupled_NNQP(i,z{j}{i},z{j}{i+1},lambda{j}{i},lambda{j}{i+1},solution, mpc_problem);
        else
            error(sprintf('Error: Unexpected input %d',aladin_options.sol_method));

        end
    end
    
    [dz{j+1}, du{j+1}, lambda{j+1}, sol_coupled{j+1}] = coupledQP(xkp_opt{j},u_opt{j}, mpc_problem,optimization_options,aladin_options);
    
    for i = 1:N+1
        z{j+1}{i} = z{j}{i} + dz{j+1}{i};
    end
    
    %% Solution result messages

    if ( aladin_options.sol_method == 0 )

        %% Implicit MPC evaluation of DecoupledQP

        for J = 1 : N
            info_decopuledQP(J) = sol_decoupled{1}{1}.problem;
        end % for J
        [info_unique_decoupledQP, info_index_unique_decoupledQP ] = unique(info_decopuledQP);

        sol_decoupled_info = [];
        for J = 1: length(info_index_unique_decoupledQP)
            sol_decoupled_info = [sol_decoupled_info, sol_decoupled{j}{ J }.info, ' | '];
        end % for J

        disp(sprintf('Iteration: %d of %d: | Coupled QP: %s | Decoupled QP: %s',j, maxiter, sol_coupled{j+1}.info, sol_decoupled_info ))
        
    elseif ( aladin_options.sol_method == 1 )

        %% Explicit MPC evaluation of DecoupledQP

        for J = 1 : 3
            info_decopuledQP(J) = solution.diagnostic{J}.problem;
        end % for J
        [info_unique_decoupledQP, info_index_unique_decoupledQP ] = unique(info_decopuledQP);

        sol_decoupled_info = [];
        for J = 1: length(info_index_unique_decoupledQP)
            sol_decoupled_info = [sol_decoupled_info,solution.diagnostic{J}.info, ' | '];
        end % for J
        
        disp(sprintf('Iteration: %d of %d: | Coupled QP: %s | Decoupled QP: %s',j, maxiter, sol_coupled{j+1}.info, sol_decoupled_info ))
    else
        error(sprintf(' Unexpected settings for ALADIN_OPTIONS.EXPLICIT: %d',aladin_options.sol_method))
    end

    j = j + 1;
    
end

sumiter = j - 1;
solution.u_opt = u_opt;
solution.xk_opt = xk_opt;
solution.lambda = lambda;
solution.sol_coupled = sol_coupled;
solution.sol_decoupled = sol_decoupled;
solution.sumiter = sumiter;
    
%% View profiles
if graph == 1

    %% Display states variables
    for j = 1:sumiter
        figure(j)
        for k = 1:1:nx
            subplot(nx, 2, k)
            hold on;box on; grid on;
            for i = 1:N
                plot([(i-1) i],[xk_opt{j}{i}(k)+xs(k) xkp_opt{j}{i}(k)+xs(k)],'x-','Color',impcolor,'LineWidth',1)
            end
            axis padded
            xlabel('prediction horizon');ylabel(['$h_' num2str(k) '$ [m]']);
            legend('ALADIN solution','Location','best')
        end
    end

    
    %% Display input variables
    for j = 1:sumiter
        u_opt_vector = cell2mat(u_opt{j});
        figure(j)
        subplot(2,2,3)
        stairs([0:N-1],[u_opt_vector(1,:)+us],'Color',impcolor,'LineWidth',linewidth,'MarkerSize',markersize)
        hold on; box on; grid on;
        axis padded;
        xlabel('prediction horizon');ylabel('$q_0$ [m$^3$/s]');
        legend('ALADIN solution','Location','southeast')
    end

end
end