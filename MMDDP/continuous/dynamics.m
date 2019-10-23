%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics
% theta=theta+dt*(theta_dot)
% theta_dot=theta_dot+dt*(sin(theta)-theta_dot+u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xnew, fx, fu]=dynamics(x,u,dt)
N=size(x,2);
n=size(x,1);

xnew=zeros(n,N);
xnew(:,1)=x(:,1);
fx(:,:,1)=[0, 1;...
        9.81*2*cos(xnew(1,1)), -0.4];
fu(:,:,1)=[0;4];

for i=1:N-1
    xnew(:,i+1)=xnew(:,i)+dt*[xnew(2,i); 9.81*2*sin(xnew(1,i))-0.4*xnew(2,i)+4*u(1,i)];
    fx(:,:,i+1)=[0, 1;...
        9.81*2*cos(xnew(1,i+1)), -0.4];
    fu(:,:,i+1)=[0;4];
end

end