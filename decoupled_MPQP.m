function [u_opt, xk_opt, xkp_opt, sol] = decoupled_MPQP(k,zK,zKP,lambdaK,lambdaKP,solution, mpc_problem)

%% Loading optimization settings

Q = mpc_problem.Q;
R = mpc_problem.R;
N = mpc_problem.N;
mpsol1 = solution.mpsol1;
mpsol2 = solution.mpsol2;
mpsol3 = solution.mpsol3;

nx = size(Q,1); 
nu = size(R,2); 

%% Defining and solving optimization problem
if (k == 1)
    theta = [zK;zKP;lambdaKP];
    [ u_out, r_found ] = my_point_location_problem( mpsol1, theta );
    u_opt = u_out(1 : nu);
    xk_opt = u_out(nu+1 : nu+nx);
    xkp_opt = u_out(nu+nx+1 : nu+nx+nx);
elseif ((k < N) & (k > 1)) 
    theta = [zK;zKP;lambdaK;lambdaKP];
    [ u_out, r_found ] = my_point_location_problem( mpsol2, theta );
    u_opt = u_out(1 : nu);
    xk_opt = u_out(nu+1 : nu+nx);
    xkp_opt = u_out(nu+nx+1 : nu+nx+nx);
elseif (k == N)
    theta = [zK;lambdaK];
    [ u_out, r_found ] = my_point_location_problem( mpsol3, theta );
    u_opt = u_out(1 : nu);
    xk_opt = u_out(nu+1 : nu+nx);
    xkp_opt = u_out(nu+nx+1 : nu+nx+nx);
else
    error(sprintf('DecoupledQP: Unexpected step k = %d', k));
end

sol = r_found;

end