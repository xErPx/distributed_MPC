function [solution] = closed_loop_simulation(mpc_problem, closed_loop_settings, optimization_options, aladin_options)

%% Display settings
linewidth = 1;
impcolor = 'blue';
mpccolor = '#808080';
% fontsize = 24;

%% Loading system specifications 
A = mpc_problem.A;
B = mpc_problem.B;
C = mpc_problem.C;
D = mpc_problem.D;
xs = [mpc_problem.x1s; mpc_problem.x2s];
us = mpc_problem.us;

%% Loading closed-loop simulation settings
Tmax = closed_loop_settings.Tmax;
N = mpc_problem.N;
Ts = closed_loop_settings.Ts;

%% Loading optimization settings
umin = mpc_problem.umin;
umax = mpc_problem.umax;
xmin = mpc_problem.xmin;
xmax = mpc_problem.xmax;

%% Simulation loop
[nx, nu] = size(B);
[ny, ~]  = size(C);

kf = Tmax/Ts;
u = zeros(nu,kf);
x = zeros(nx,kf);
y = zeros(ny,kf);
u_mpc = zeros(nu,kf);
x_mpc = zeros(nx,kf);
y_mpc = zeros(ny,kf);
x(:,1) = mpc_problem.x0;
x_mpc(:,1) = mpc_problem.x0;

for k = 1:1:kf

    disp(sprintf('*** Closed-loop simulation step: %d of %d: | (ALADIN) ***',k, kf))
    solution = aladin_solver(mpc_problem, optimization_options, aladin_options);
    
    sumiter = solution.sumiter;
    u(:,k) = solution.u_opt{end}{1}; %solution.u_opt{iteration}{step_of_prediction_horizon}

    x(:, k+1) = A*x(:, k) + B*u(:, k);
    y(:, k) = C*x(:, k) + D*u(:, k);

    mpc_problem.x0 = x(:, k+1);

end

mpc_problem.x0 = x_mpc(:, 1);
for k = 1:1:kf

    disp(sprintf('*** Closed-loop simulation step: %d of %d: | (MPC)***',k, kf))

    solution = ModelPredictiveControl(mpc_problem,optimization_options,solution);
    
    u_mpc(:,k) = solution.u_mpc{1}; 

    x_mpc(:, k+1) = A*x_mpc(:, k) + B*u_mpc(:, k);
    y_mpc(:, k) = C*x_mpc(:, k) + D*u_mpc(:, k);

    mpc_problem.x0 = x_mpc(:, k+1);

end

x(:,k+1) = [];
x_mpc(:,k+1) = [];

%% Views

path = [cd, '\Figures\spojite'];
time = linspace(0, kf*Ts, kf);
figure
for k = 1:1:nx
    subplot(nx, 1, k)
    hold on
    axis padded;

    stairs(time, x(k, :) + xs(k), ...
        'linewidth', linewidth, 'Color', impcolor);
    stairs(time, x_mpc(k, :) + xs(k), ...
        'LineStyle', '--', 'Color', mpccolor, 'linewidth', linewidth);

    xlabel('$t$ [s]')
    ylabel(['$h_' num2str(k) '$ [m]'])
    set(gcf, 'name', 'states')
%     name = sprintf('CL_states_15');
%     set(gcf,'PaperOrientation','landscape');
%     set(gcf,'PaperUnits','normalized')
%     set(gcf,'PaperPosition',[0 0 1 1])
%     set(gca,'fontsize',fontsize)
%     saveas(gcf,fullfile(path, [name,'.pdf']));
end

figure
for k = 1:1:nu
    subplot(nu, 1, k)
    hold on
    axis padded;
    
    stairs(time, u(k, :) + us(k), ...
        'linewidth', linewidth, 'Color', impcolor);
    stairs(time, u_mpc(k, :) + us(k), ...
        'LineStyle', '--', 'Color', mpccolor, 'linewidth', linewidth);
    xlabel('$t$ [s]')
    ylabel(['$q_0$ [m$^3$/s]'])
    set(gcf, 'name', 'inputs')
%     name = sprintf('CL_input_15');
%     set(gcf,'PaperOrientation','landscape');
%     set(gcf,'PaperUnits','normalized')
%     set(gcf,'PaperPosition',[0 0 1 1])
%     set(gca,'fontsize',fontsize)
%     saveas(gcf,fullfile(path,[name,'.pdf']));
end

figure
hold on
axis padded;

stairs(time, y + xs(2), ...
    'linewidth', linewidth, 'Color', impcolor);
stairs(time, y_mpc + xs(2), ...
    'LineStyle', '--', 'Color', mpccolor, 'linewidth', linewidth);
yline(xs(2), 'linestyle', '--', 'linewidth', linewidth, 'Color', 'k');
xlabel('$t$ [s]')
ylabel(['$h_2$ [m]'])
set(gcf, 'name', 'outputs')
% name = sprintf('CL_output_15');
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperUnits','normalized')
% set(gcf,'PaperPosition',[0 0 1 1])
% set(gca,'fontsize',fontsize)
% saveas(gcf,fullfile(path,[name,'.pdf']));
end