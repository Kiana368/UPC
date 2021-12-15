%function [vr,vt,time,loss,radius]=escape()

G=6.67*10^(-11); %gravitational 
v_radius=0;
M=5.965*10^(24);
m0=5000;       %initial mass
r(1)=6.65*10^6; %initial radius
w_Theta=7735/(r(1));
F=0.4;     
A=1.02*10^(-5); %rate of mass loss
t=0;          
dt=0.01;          %time interval
theta(1)=0;       %initial angle (assumed as zero)
index=1;          
while(1)
   [a_radius_tem1,b_theta_tem1]=acceleration(v_radius,w_Theta,t,r(index),m0);
   v_radius_tem= v_radius+a_radius_tem1*dt;
   v_theta_tem= w_Theta+b_theta_tem1*dt;
   [a_radius_tem2,b_theta_tem2]=acceleration(v_radius_tem,v_theta_tem,t,r(index),m0-A*dt);
   v_radius=v_radius+1/2*(a_radius_tem1+a_radius_tem2)*dt;
   w_Theta=w_Theta+1/2*(b_theta_tem1+b_theta_tem2)*dt;
   index=index+1;
   theta(index)=theta(index-1)+w_Theta*dt;
   r(index)=r(index-1)+v_radius*dt;
   m0=m0-A*dt;
   t=t+dt;
   e=energy(v_radius,w_Theta,r(index));
   if (e>=0)
       break;
   end
   
end
polar(theta,r);
grid on;
%fprintf("Time for escape: %f s\n",t);
%fprintf("loss of mass %f kg",A*t);
time=t;
loss=A*t;
radius=r(index);
%disp(r);
%end           
