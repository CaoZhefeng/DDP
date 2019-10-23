%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics
% theta=theta+dt*(theta_dot)
% theta_dot=theta_dot+dt*(sin(theta)-theta_dot+u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xnew, dx, fx, fu, fv]=dynamics(x,u,v,index,dt)
N=size(x,2);
n=size(x,1);
m=size(u,1);

if index==0
    xnew=zeros(n,N);
    xnew(:,1)=x(:,1);
    fx(:,:,1)=[0, 1;...
            9.81*2*cos(xnew(1,1)), -0.4]*dt+eye(n);
    fu(:,:,1)=[0;4]*dt;
    fv(:,:,1)=[0;4]*dt;

    for i=1:N-1
        xnew(:,i+1)=xnew(:,i)+dt*[xnew(2,i); 9.81*2*sin(xnew(1,i))-0.4*xnew(2,i)+4*u(1,i)+4*v(1,i)];
        fx(:,:,i+1)=[0, 1;...
            9.81*2*cos(xnew(1,i+1)), -0.4]*dt+eye(n);
        fu(:,:,i+1)=[0;4]*dt;
        fv(:,:,i+1)=[0;4]*dt;
    end
    dx=xnew-x;
else
    xnew=x(:,index)+dt*[x(2,index); 9.81*2*sin(x(1,index))-0.4*x(2,index)+4*u(1,index)+4*v(1,index)];
    fx=[0, 1;...
        9.81*2*cos(x(1,index)), -0.4]*dt+eye(n);
    fu=[0;4]*dt;
    fv=[0;4]*dt;
    dx=xnew-x(:,index+1);
end

end