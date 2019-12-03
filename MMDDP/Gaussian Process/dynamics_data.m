function dxdt=dynamics_data(t,x,u_pp)
dxdt=zeros(4,1);
u = ppval(u_pp, t);

m=1;
l=0.5;
g=9.81;

u = u+0.01*u.*randn(2,1);

dxdt(1)=x(3);
dxdt(2)=x(4);
dxdt(3)=((u(1)-u(2)*cos(x(1)-x(2)))/(m*l^2)-sin(x(1)-x(2))*(x(3)^2*cos(x(1)-x(2))+x(4)^2)+g/l*(2*sin(x(1))-sin(x(2))*cos(x(1)-x(2))))/(2-cos(x(1)-x(2))^2);
dxdt(4)=((2*u(2)-u(1)*cos(x(1)-x(2)))/(m*l^2)+sin(x(1)-x(2))*(2*x(3)^2+x(4)^2*cos(x(1)-x(2)))+2*g/l*(sin(x(2))-sin(x(1))*cos(x(1)-x(2))))/(2-cos(x(1)-x(2))^2);
end