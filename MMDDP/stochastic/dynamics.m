%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics
% theta=theta+dt*(theta_dot)
% theta_dot=theta_dot+dt*(sin(theta)-theta_dot+u+v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdt=dynamics(t,x,u_pp,v_pp)
u = ppval(u_pp, t);
v = ppval(v_pp, t);

dxdt=[x(2);9.81*2*sin(x(1))-0.4*x(2)+4*u+4*v];
end
