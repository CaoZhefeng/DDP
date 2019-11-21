% vectorized dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics
% theta=theta+dt*(theta_dot)
% theta_dot=theta_dot+dt*(sin(theta)-theta_dot+u+v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dxdt is n*d matrix n is the number of data d is the dimension of state
% u v is both n*k k is the dimension of control
% consider torque stochastic

function dxdt=duopend(x,u,v)

m=1;
l=0.5;
g=9.81;
n = size(x,1);

u = 0.95*u+0.1*u.*rand(n,2);
v = 0.95*v+0.1*v.*rand(n,2);

dxdt(:,1)=x(:,3);
dxdt(:,2)=x(:,4);
dxdt(:,3)=(((u(:,1)+v(:,1))-(u(:,2)+v(:,2)).*cos(x(:,1)-x(:,2)))./(m*l^2)-sin(x(:,1)-x(:,2)).*(x(:,3).^2.*cos(x(:,1)-x(:,2))+x(:,4).^2)+g/l*(2*sin(x(:,1))-sin(x(:,2)).*cos(x(:,1)-x(:,2))))./(2-cos(x(:,1)-x(:,2)).^2);
dxdt(:,4)=((2*(u(:,2)+v(:,2))-(u(:,1)+v(:,1)).*cos(x(:,1)-x(:,2)))./(m*l^2)+sin(x(:,1)-x(:,2)).*(2*x(3).^2+x(4).^2.*cos(x(:,1)-x(:,2)))+2*g/l*(sin(x(:,2))-sin(x(:,1)).*cos(x(:,1)-x(:,2))))./(2-cos(x(:,1)-x(:,2)).^2);
end
