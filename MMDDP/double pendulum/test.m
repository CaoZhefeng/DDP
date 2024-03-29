syms m l g
x = sym('x',[4,1]);
u = sym('u',[2,1]);
f(1)=x(3);
f(2)=x(4);
f(3)=((u(1)-u(2)*cos(x(1)-x(2)))/(m*l^2)-sin(x(1)-x(2))*(x(3)^2*cos(x(1)-x(2))+x(4)^2)+g/l*(2*sin(x(1))-sin(x(2))*cos(x(1)-x(2))))/(2-cos(x(1)-x(2))^2);
f(4)=((2*u(2)-u(1)*cos(x(1)-x(2)))/(m*l^2)+sin(x(1)-x(2))*(2*x(3)^2+x(4)^2*cos(x(1)-x(2)))+2*g/l*(sin(x(2))-sin(x(1))*cos(x(1)-x(2))))/(2-cos(x(1)-x(2))^2);

fx=[0 0 1 0;...
    0 0 0 1;...
    diff(f(3),x(1)), diff(f(3),x(2)), diff(f(3),x(3)), diff(f(3),x(4));...
    diff(f(4),x(1)), diff(f(4),x(2)), diff(f(4),x(3)), diff(f(4),x(4))];

fu=[0 0;...
    0 0;...
    diff(f(3),u(1)), diff(f(3),u(2));...
    diff(f(4),u(1)), diff(f(4),u(2))];