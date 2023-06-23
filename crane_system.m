function dXdt = crane_system(t,X,p)

    % Extraction of states
    w = X(1:p.N);
    w_t = X(p.N+1:2*p.N);
    phi = X(2*p.N+1);
    phi_t = X(2*p.N+2);
    z = X(2*p.N+3);
    z_t = X(2*p.N+4);

    %% Controler 
    epsilon = p.z_soll-z;
    if p.controller_bool
        % Controller
        F_fr = calculate_friction_force(z_t,p);
        F_t = controller_output(z_t,F_fr,epsilon,p);
    else
        % No controller
        F_fr = 0;
        F_t = 0;
    end

    

    %% ODE for phi
    phi_tt = -1/(p.l * p.m_k) * (p.g*p.m_1*phi + F_t - F_fr);
   
    %% ODE for z
    b_b = p.rho*p.delta_x/(p.rho*p.delta_x + 2*p.m_ges);                                                            

    z_tt = b_b*p.EI/(p.rho*p.delta_x^4)*(2*w(p.N-2) - 4*w(p.N-1) + 2*w(p.N)) + ...                                  
           b_b*p.c/p.rho*w_t(p.N) + ...                                                                             
           p.g*p.m_l/p.m_k*phi + ...
           (1/p.m_k + 2/(2*p.m_t+p.rho*p.delta_x))*(F_t - F_fr);                                                    
    
    %% ODEs for w_tt
    w_tt(1) = 1/p.rho*(p.EI * ...
              ((7*w(1)-4*w(2)+w(3))/(p.delta_x^4)) - p.c*w_t(1));
    w_tt(2) = 1/p.rho*(p.EI * ...
              ((-4*w(1)+6*w(2)-4*w(3)+w(4))/(p.delta_x^4)) - p.c*w_t(2));
    for i=3:p.N-2
        w_tt(i) = 1/p.rho * (p.EI * ...
                  ((w(i-2) - 4*w(i-1) + 6*w(i) - 4*w(i+1) + w(i+2))/ ...
                  (p.delta_x^4)) - p.c*w_t(i));
    end

    w_tt(p.N-1) = 1/p.rho * (p.EI * ...
                  ((w(p.N-3)-4*w(p.N-2)+5*w(p.N-1)-2*w(p.N))/(p.delta_x^4)) - p.c*w_t(p.N-1));

    w_tt(p.N) = -b_b*p.EI/p.rho/(p.delta_x^4)*(2*w(p.N) - 4*w(p.N-1) + 2*w(p.N-2)) ...
                -b_b*p.c/p.rho*w_t(p.N) + ...
                -b_b*2/p.rho/p.delta_x*(F_t - F_fr);

    %% Return
    dXdt = [w_t',w_tt,phi_t,phi_tt,z_t,z_tt]';

    [t F_t epsilon w_tt(p.N) w(end) p.z_soll]
end