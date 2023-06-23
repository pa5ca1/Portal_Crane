clc; clear; close all;

%% Parameter
p.g = 9.81; % [m/s^2] gravitation
p.rho = 6124; % [kg/m^3] 
p.L = 0.5; % [m]
p.m_t = 1.46; % [kg]
p.m_k = .1; % [kg]
p.m_l = .1; % [kg]
p.l = .67; % [m]
p.EI = 4340; % []
p.c = .05; % []

% Self-defined parameters
p.m_1 = p.m_k + p.m_l;
p.m_ges = p.m_k + p.m_l + p.m_t;
% Number of support points
p.N = 20;
p.delta_x = p.L/p.N;
% Controler
p.epsilon = 0.1;
p.k = 1;
p.k_p = 1;
p.d = 1;
p.z_soll = 1;


%% ODE solver
% Initical conditions
x0 = zeros(1,2*p.N+4);
% Diviosn by w_t(N)
x0(2*p.N) = 0.1;
%x0(2*p.N+3) = 2;
% Simulation time
tspan = [0 50];
% Solving
options = [];
[t,X] = ode113(@(t,X) crane_system(t,X,p),tspan,x0,options);

figure(1)
plot(t,X(:,2*p.N+1),'+-','Color','red')
hold on;
plot(t,X(:,2*p.N+2),'+-','Color','blue')
hold off;
legend('\varphi','\varphi_t');

figure(2)
plot(t,X(:,2*p.N+3),'+-','Color','red')
hold on;
plot(t,X(:,2*p.N+4),'+-','Color','blue')
legend();


function dXdt = crane_system(~,X,p)
    % Extraction of states
    % Are these really the states?
    w = X(1:p.N);
    w_t = X(p.N+1:2*p.N);
    phi = X(2*p.N+1);
    phi_t = X(2*p.N+2); % unused?
    z = X(2*p.N+3);
    z_t = X(2*p.N+4);

    %% Controler ???
    % Calculation epsilon
    epsilon = z-p.z_soll;
    F_t = (-p.d*z_t + p.k_p * epsilon);
    %F_t = 0.01;
    
    % Ansatz 1: F_fr = 0
    % Ansatz 2: F_fr = p.mu * z_t
    F_fr = 0;

    %% ODE for phi
    %phi_tt = -1/(p.l * p.m_k) * (p.g*p.m_1*phi + F_t * (p.mu-1));
    phi_tt = -1/(p.l * p.m_k) * (p.g*p.m_1*phi + F_fr - F_t);
   
    %% ODE for z
    b_b = p.rho*p.delta_x/(p.rho*p.delta_x + 2*p.m_t); % can be a parameter in structure p.b_b
    z_tt = b_b * p.EI/p.rho*p.delta_x^4 * (2*w(p.N-2) - 4*w(p.N-1) + 2*w(p.N)) + ...
           b_b * p.c/(p.rho * w_t(p.N)) + ...
           p.g*p.m_l/p.m_k * (1/p.m_k + 2/(2*p.m_k+p.rho*p.delta_x)) * (-F_t + F_fr);
    
    %% ODEs for w_tt
    % i=1: We have a problem with w(i-2)=w(-1) and w(i-1)=w(0)
    % Solution: Plug in w(-1)=w(1) and w(0)=0
    w_tt(1) = 1/p.rho * (p.EI * ...
              ((w(1)+6*w(1)-4*w(2)+w(3))/(p.delta_x^4)) - p.c*w_t(1));
    % i=2: We still have a problem with w(i-1)=w(0)
    % Solution: Plug in w(0)=0
    w_tt(2) = 1/p.rho * (p.EI * ...
              ((4*w(1)+6*w(2)-4*w(3)+w(4))/(p.delta_x^4)) - p.c*w_t(2));
    for i=3:p.N-2
        w_tt(i) = 1/p.rho * (p.EI * ...
                  ((w(i-2) - 4*w(i-1) + 6*w(i) - 4*w(i+1) + w(i+2))/ ...
                  (p.delta_x^4)) - p.c*w_t(i));
    end
    % i=N-1: We have a problem with w(i+2)=w(N+1)
    % Solution: Plug in  w(N+1)=2*w(N)-w(N-1)
    w_Np1 = 2*w(p.N) - w(p.N-1); % w(N+1)
    w_tt(p.N-1) = 1/p.rho * (p.EI * ...
                  ((w(p.N-3)-4*w(p.N-2)+6*w(p.N-1)-4*w(p.N)+ w_Np1 )/(p.delta_x^4)) - p.c*w_t(i));
    % i=N: See project work
    w_tt(p.N) = -p.EI/(p.rho*p.delta_x^4 + 2*p.delta_x^3*p.m_ges) * ... 
                 ( 2*w(p.N) - 4*w(p.N-1) + 2*w(p.N-2) + ...
                 2*p.delta_x^3/p.EI * (p.m_l*p.l*phi_tt + p.m_1*z_tt) - ...
                 (p.c*p.delta_x^4/(p.rho*p.EI)) * w_t(p.N) );

    %% Return
    % ??? 2*p.N+4 states but only p.N+2 returns ???
    dXdt = [w_t',w_tt,phi_t,phi_tt,z_t,z_tt]';
end












