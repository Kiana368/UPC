function [a_radius,b_theta]=acceleration(v_radius,w_theta,t,r,m0)         % v_radius :the increasing rate of radius,  v_theta:angular speed
   G=6.67*10^(-11);         %gravitational constant
   M=5.965*10^(24);         %mass of earth
   F=0.4;                   %constant force 
   A=1.02*10^(-5);          %rate of mass loss
   a_radius=-G*M/(r^2)+v_radius/(sqrt((r*w_theta)^2+v_radius^2))*F/(m0-A*t)+w_theta^2*r;
   b_theta=(r*w_theta/(sqrt((r*w_theta)^2+v_radius^2))*F/(m0-A*t)-2*v_radius*w_theta)/r;
   
end