clc; clear all; close all; yalmip clear;

%% Process inicialization

%Constants
F1 = 0.8;
F2 = 2.0;
k11 = 1.0;
k22 = 1.3;
q0s = 1;
Ts = 1.5;

%Steady states
h1s = (q0s/k11)^2;
h2s = (2*q0s/k22)^2;

%State-space model
K1 = k11/(2*sqrt(h1s));
K2 = k22/(2*sqrt(h2s));

Ac = [-K1/F1    0;
      K1/F2 -K2/F2];
Bc = [1/F1; 1/F2];
Cc = [0 1];
Dc = [0];
sysc= ss(Ac,Bc,Cc,Dc);

%Discrete time state-space model
sysd = c2d(sysc,Ts);

A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

%Process constrains
x0 = [-0.3; 2]; %zaciatocne podmienky x(k=0)
umin = -0.1; %ohranicenie riadenia
umax = 0.1;
xmin = [-0.5; -0.5]; %ohranicenie stavov
xmax = [0.5; 2.5];

%MPC tuning parameters
Q = 15*diag(1./([h1s; h2s].^2));
R = 10*diag(1./(q0s.^2));
N = 20; %prediction horizon

%% Setting algorithm parameters

solver = 'gurobi';
solver_verbose = 0; % 0 -> minimal display; 1 -> modest display level

%ALADIN tuning parameters
rho = 5; %penalty parameter
maxiter = 15; %maximum number of iterations
tolerance = 0.01;
init_var = 'dlqr'; %initialization of omega (scaling matrix)
% init_var = 'zeros';

%Display settings
% graph = 0; %solution without graphs
graph = 1; %solution with graphs

%Solution method
% sol_method = 1; %explicit MPC
sol_method = 0; %implicit MPC

%% Simulation settings
Tmax = 30; %time of simulation
% closed_loop = 1; %run closed-loop simulation
closed_loop = 0; %run prediction of trajectory

%% Creating structures for solver function

mpc_problem.A = A;
mpc_problem.B = B;
mpc_problem.C = C;
mpc_problem.D = D;
mpc_problem.umin = umin;
mpc_problem.umax = umax;
mpc_problem.xmin = xmin;
mpc_problem.xmax = xmax;
mpc_problem.Q = Q;
mpc_problem.R = R;
mpc_problem.x0 = x0;
mpc_problem.N = N;
mpc_problem.x1s = h1s;
mpc_problem.x2s = h2s;
mpc_problem.us= q0s;

closed_loop_settings.Tmax = Tmax;
closed_loop_settings.Ts = Ts;

optimization_options.yalmip.solver = solver;
optimization_options.yalmip.verbose = 0;
optimization_options.graph = graph;

aladin_options.rho = rho;
aladin_options.maxiter = maxiter;
aladin_options.init_var = init_var;
aladin_options.tolerance = tolerance;
aladin_options.sol_method = sol_method;

%% Calling the solver function

if closed_loop == 0
    solution = aladin_solver(mpc_problem, optimization_options, aladin_options);
    solution = ModelPredictiveControl(mpc_problem,optimization_options,solution);
elseif closed_loop == 1
    solution = closed_loop_simulation(mpc_problem, closed_loop_settings, optimization_options, aladin_options);
end