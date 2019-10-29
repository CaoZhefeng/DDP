%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics
% theta=theta+dt*(theta_dot)
% theta_dot=theta_dot+dt*(sin(theta)-theta_dot+u+v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdt=dynamics(t,x,ut,u,vt,v)
dxdt=zeros(2,1);
u = interp1(ut,u,t);
v = interp1(vt,v,t);


dxdt(1)=x(2);
dxdt(2)=9.81*2*sin(x(1))-0.4*x(2)+4*u+4*v;
end
