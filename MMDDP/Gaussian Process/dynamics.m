%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics
% theta=theta+dt*(theta_dot)
% theta_dot=theta_dot+dt*(sin(theta)-theta_dot+u+v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdt=dynamics(t,x,u_pp,v_pp,gprMd1,gprMd2,gprMd3,gprMd4)
dxdt=zeros(4,1);
u = ppval(u_pp, t);
v = ppval(v_pp, t);
m=1;
l=0.5;
g=9.81;

dxdt(1)=predict(gprMd1,[x',u']);
dxdt(2)=predict(gprMd2,[x',u']);
dxdt(3)=predict(gprMd3,[x',u'])+v(1);
dxdt(4)=predict(gprMd4,[x',u'])+v(2);

A(1)=x(3);
A(2)=x(4);
A(3)=(((u(1)+v(1))-(u(2)+v(2))*cos(x(1)-x(2)))/(m*l^2)-sin(x(1)-x(2))*(x(3)^2*cos(x(1)-x(2))+x(4)^2)+g/l*(2*sin(x(1))-sin(x(2))*cos(x(1)-x(2))))/(2-cos(x(1)-x(2))^2);
A(4)=((2*(u(2)+v(2))-(u(1)+v(1))*cos(x(1)-x(2)))/(m*l^2)+sin(x(1)-x(2))*(2*x(3)^2+x(4)^2*cos(x(1)-x(2)))+2*g/l*(sin(x(2))-sin(x(1))*cos(x(1)-x(2))))/(2-cos(x(1)-x(2))^2);
dxdt(1)-A(1)
dxdt(2)-A(2)
dxdt(3)-A(3)
dxdt(4)-A(4)
end