clc; clear; close all;

%% Parameter
p.g = 9.81; % [m/s^2] gravitation
p.rho = 6124; % [kg/m^3] 
p.L = 0.5; % [m]
p.m_t = 1.46; % [kg]
p.m_k = .1; % [kg]
p.m_l = .1; % [kg]
p.l = .67; % [m]
p.EI = 7e9; % []
p.c = .05; % []
p.mu = 0.5; % Friction coefficient


% Self-defined parameters
p.m_1 = p.m_k + p.m_l;
p.m_ges = p.m_k + p.m_l + p.m_t;
% Number of support points
p.N = 20;
p.delta_x = p.L/(p.N-1);

%% Controler
p.controller_bool = 1;

p.k = 0.01;
p.k_p = 0.01;
p.d = 0.01;

p.z_soll = 1;


%% ODE solver
% Initical conditions
x0 = zeros(1,2*p.N+4);
% Simulation time
t_end = 300;
tspan = 0:0.01:t_end;
% Solving
options = odeset('RelTol',1e-5,'AbsTol',1e-7);
%[t,X] = ode23s(@(t,X) crane_system(t,X,p),tspan,x0,options);
[t,X] = ode15s(@(t,X) crane_system(t,X,p),tspan,x0,options);


w = X(:,1:p.N);
w_t = X(:,p.N+1:2*p.N);
phi = X(:,2*p.N+1);
phi_t = X(:,2*p.N+2);
z = X(:,2*p.N+3);
z_t = X(:,2*p.N+4);

figure(1)
plot(t,phi_t,'+-','Color','blue','DisplayName','Angular velocity (phi_t)')
hold on;
plot(t,phi,'+-','Color','red','DisplayName','Angle (phi)')
hold off;
ylabel('Angle and Angular velocity')
xlabel('Time t')
%legend('varphi','varphi_t');
legend();
title('Angle and angular velocity')

figure(2)
plot(t,z,'+-','Color','red','DisplayName','Position (z)')
hold on;
plot(t,z_t,'+-','Color','blue','DisplayName','Velocity (z_t)')
xlabel('Time t')
ylabel('Position and velocity')
legend();
title('Position and velocity of trolley')

figure(3)
[t_grid,x_grid] = meshgrid(tspan,0:p.delta_x:p.L);
surf(t_grid',w,x_grid',w,'EdgeColor','none');
title('Beam deflection')
xlabel('Time t')
ylabel('Deflection w')
zlabel('Height x')
colormap('jet');
clim([min(min(w)) max(max(w))])
colorbar();
view([25.1450,42.0175])

figure(4)
[t_grid,x_grid] = meshgrid(tspan,0:p.delta_x:p.L);
surf(t_grid',w_t,x_grid',w_t,'EdgeColor','none');
title('Beam deflection velocity')
xlabel('Time t')
ylabel('Deflection w')
zlabel('Height x')
colormap('jet');
clim([min(min(w_t)) max(max(w_t))])
colorbar();
view([25.1450,42.0175])

figure(5)
epsilon = p.z_soll-z;
F_fr = calculate_friction_force(z_t,p);
F_t = controller_output(z_t,F_fr,epsilon,p);
plot(t,F_t)
title('Contoller output u')














