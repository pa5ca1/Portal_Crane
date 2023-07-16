clc; clear; close all;


%% Make everything latex
% Source: https://de.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
set(groot,'defaultAxesFontSize',24)

%% Parameter
p.g = 9.81;     % [m/s^2] gravitation
p.rho = 6124;   % [kg/m^3] density 
p.L = 0.5;      % [m] length of crane rope
p.m_t = 1.46;   % [kg] mass of bridge
p.m_k = .1;     % [kg] mass of gantry
p.m_l = .1;     % [kg] mass of container
p.l = .67;      % [m] length between mass and gantry
p.EI = 7e9;     % [Nm^2] Young's modulus
p.c = .05;      % [1] dissipation constant
p.mu = 0.5;     % [1] friction constant
% Self-defined parameters
p.m_1 = p.m_k + p.m_l;
p.m_ges = p.m_k + p.m_l + p.m_t;
% Number of support points
p.N = 20;
p.delta_x = p.L/(p.N-1);

%% Controler
% Choose if we want to have a controller (1) or not (0)
p.controller_bool = 1;
% Values for controller parameters
p.k = 0.01;
p.k_p = 0.01;
p.d = 0.01;
% reference value for position
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

%% Plots

figure(1)
plot(t,phi_t,'+-','Color','blue','DisplayName','Angular velocity ($\varphi_t$)')
hold on;
plot(t,phi,'+-','Color','red','DisplayName','Angle ($\varphi$)')
hold off;
ylabel('Angle and Angular velocity')
xlabel('Time t')
%legend('varphi','varphi_t');
legend();
%title('Angle and angular velocity')

figure(2)
plot(t,z,'+-','Color','red','DisplayName','Position ($z$)')
hold on;
plot(t,z_t,'+-','Color','blue','DisplayName','Velocity ($z_t$)')
xlabel('Time $t$')
ylabel('Position $z$ and velocity $z_t$')
legend();
%title('Position and velocity of trolley')

fig = figure(3);
fig.Position = [100 100 950 950];
[t_grid,x_grid] = meshgrid(tspan,0:p.delta_x:p.L);
surf(t_grid',w,x_grid',w,'EdgeColor','none');
%title('Beam deflection $w(t,x)$')
xlabel('Time $t$')
ylabel('Deflection $w$')
zlabel('Height $x$')
colormap('jet');
try
    clim([min(min(w)) max(max(w))])
catch
    
end
colorbar('NorthOutside');
view([25.1450,42.0175])
orient(fig,'landscape')
set(gcf,'PaperPositionMode','auto')
print('beam_deflection','-dpdf','-r0','-fillpage')


fig = figure(4);
fig.Position = [100 100 950 950];
[t_grid,x_grid] = meshgrid(tspan,0:p.delta_x:p.L);
surf(t_grid',w_t,x_grid',w_t,'EdgeColor','none');
%title('Beam deflection velocity $w_t(t,x)$')
xlabel('Time $t$')
ylabel('Deflection velocity $w_t$')
zlabel('Height $x$')
colormap('jet');
try
clim([min(min(w_t)) max(max(w_t))])
catch

end
colorbar("NorthOutside");
view([25.1450,42.0175])
orient(fig,'landscape')
set(gcf,'PaperPositionMode','auto')
print('beam_velocity','-dpdf','-r0','-fillpage')

figure(5)
epsilon = p.z_soll-z;
F_fr = calculate_friction_force(z_t,p);
F_t = controller_output(z_t,F_fr,epsilon,p);
plot(t,F_t,'+-','Color','blue')
xlabel('Time $t$')
ylabel('Contoller output $u$')




figure(6)
t0=9000;
t1=11500;
plot(t(t0:t1),phi_t(t0:t1),'+-','Color','blue','DisplayName','Angular velocity ($\varphi_t$)')
hold on;
plot(t(t0:t1),phi(t0:t1),'+-','Color','red','DisplayName','Angle ($\varphi$)')
hold off;
ylabel('Angle and Angular velocity')
xlabel('Time t')
xlim([t(t0),t(t1)])
legend()
%title('Position and velocity of trolley')









