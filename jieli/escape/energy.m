function e=energy(v_radius,w_theta,r)
    G=6.67*10^(-11);
    M=5.965*10^(24);
    e=1/2*(v_radius^2+(w_theta*r)^2)-G*M/r;
    disp(e);
end