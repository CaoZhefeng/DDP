%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics
% theta=theta+dt*(theta_dot)
% theta_dot=theta_dot+dt*(sin(theta)-theta_dot+u+v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdt=dynamics(t,x,u_pp,v_pp)
dxdt=zeros(2,1);
u = ppval(u_pp, t);
v = ppval(v_pp, t);


dxdt(1)=x(2);
dxdt(2)=9.81*2*sin(x(1))-0.4*x(2)+4*u+4*v;
end
