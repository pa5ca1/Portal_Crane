function F_t = controller_output(z_t,F_fr,epsilon,p)
    F_t = -p.d*z_t + F_fr + p.k_p * epsilon;
end