function dVdt = odeback(t,V,u_pp,v_pp,x_pp,Q,Ru,Rv,alpha)

V_x = V(1:4);
V_xx=reshape(V(5:end), 4, 4);
u = ppval(u_pp, t);
v = ppval(v_pp, t);
x = ppval(x_pp, t);
m=1;
l=0.5;
g=9.81;

fx=[                                                                                                                                                                                                                                                                                                                                                             0,                                                                                                                                                                                                                                                                                                                                                   0,                                                     1,                                                      0;...
                                                                                                                                                                                                                                                                                                                                                             0,                                                                                                                                                                                                                                                                                                                                                   0,                                                     0,                                                      1;...
                    - (x(3)^2*sin(x(1) - x(2))^2 - cos(x(1) - x(2))*(cos(x(1) - x(2))*x(3)^2 + x(4)^2) + (g*(2*cos(x(1)) + sin(x(1) - x(2))*sin(x(2))))/l + ((u(2)+v(2))*sin(x(1) - x(2)))/(l^2*m))/(cos(x(1) - x(2))^2 - 2) - (2*cos(x(1) - x(2))*sin(x(1) - x(2))*((g*(2*sin(x(1)) - cos(x(1) - x(2))*sin(x(2))))/l - sin(x(1) - x(2))*(cos(x(1) - x(2))*x(3)^2 + x(4)^2) + ((u(1)+v(1)) - (u(2)+v(2))*cos(x(1) - x(2)))/(l^2*m)))/(cos(x(1) - x(2))^2 - 2)^2, (x(3)^2*sin(x(1) - x(2))^2 - cos(x(1) - x(2))*(cos(x(1) - x(2))*x(3)^2 + x(4)^2) + (g*(cos(x(1) - x(2))*cos(x(2)) + sin(x(1) - x(2))*sin(x(2))))/l + ((u(2)+v(2))*sin(x(1) - x(2)))/(l^2*m))/(cos(x(1) - x(2))^2 - 2) + (2*cos(x(1) - x(2))*sin(x(1) - x(2))*((g*(2*sin(x(1)) - cos(x(1) - x(2))*sin(x(2))))/l - sin(x(1) - x(2))*(cos(x(1) - x(2))*x(3)^2 + x(4)^2) + ((u(1)+v(1)) - (u(2)+v(2))*cos(x(1) - x(2)))/(l^2*m)))/(cos(x(1) - x(2))^2 - 2)^2, (2*x(3)*cos(x(1) - x(2))*sin(x(1) - x(2)))/(cos(x(1) - x(2))^2 - 2),               (2*x(4)*sin(x(1) - x(2)))/(cos(x(1) - x(2))^2 - 2);...
 - (cos(x(1) - x(2))*(2*x(3)^2 + cos(x(1) - x(2))*x(4)^2) - x(4)^2*sin(x(1) - x(2))^2 - (2*g*(cos(x(1) - x(2))*cos(x(1)) - sin(x(1) - x(2))*sin(x(1))))/l + ((u(1)+v(1))*sin(x(1) - x(2)))/(l^2*m))/(cos(x(1) - x(2))^2 - 2) - (2*cos(x(1) - x(2))*sin(x(1) - x(2))*(sin(x(1) - x(2))*(2*x(3)^2 + cos(x(1) - x(2))*x(4)^2) + (2*(u(2)+v(2)) - (u(1)+v(1))*cos(x(1) - x(2)))/(l^2*m) + (2*g*(sin(x(2)) - cos(x(1) - x(2))*sin(x(1))))/l))/(cos(x(1) - x(2))^2 - 2)^2,      (cos(x(1) - x(2))*(2*x(3)^2 + cos(x(1) - x(2))*x(4)^2) - x(4)^2*sin(x(1) - x(2))^2 - (2*g*(cos(x(2)) - sin(x(1) - x(2))*sin(x(1))))/l + ((u(1)+v(1))*sin(x(1) - x(2)))/(l^2*m))/(cos(x(1) - x(2))^2 - 2) + (2*cos(x(1) - x(2))*sin(x(1) - x(2))*(sin(x(1) - x(2))*(2*x(3)^2 + cos(x(1) - x(2))*x(4)^2) + (2*(u(2)+v(2)) - (u(1)+v(1))*cos(x(1) - x(2)))/(l^2*m) + (2*g*(sin(x(2)) - cos(x(1) - x(2))*sin(x(1))))/l))/(cos(x(1) - x(2))^2 - 2)^2,             -(4*x(3)*sin(x(1) - x(2)))/(cos(x(1) - x(2))^2 - 2), -(2*x(4)*cos(x(1) - x(2))*sin(x(1) - x(2)))/(cos(x(1) - x(2))^2 - 2)];

fu=[                                         0,                                         0;...
                                         0,                                         0;...
           -1/(l^2*m*(cos(x(1) - x(2))^2 - 2)), cos(x(1) - x(2))/(l^2*m*(cos(x(1) - x(2))^2 - 2));...
 cos(x(1) - x(2))/(l^2*m*(cos(x(1) - x(2))^2 - 2)),           -2/(l^2*m*(cos(x(1) - x(2))^2 - 2))];
fv=[                                         0,                                         0;...
                                         0,                                         0;...
           -1/(l^2*m*(cos(x(1) - x(2))^2 - 2)), cos(x(1) - x(2))/(l^2*m*(cos(x(1) - x(2))^2 - 2));...
 cos(x(1) - x(2))/(l^2*m*(cos(x(1) - x(2))^2 - 2)),           -2/(l^2*m*(cos(x(1) - x(2))^2 - 2))];

F=[0;0;alpha*u(1)^2;alpha*u(2)^2];
F_x=zeros(4,4);
F_u=[0, 0;0, 0;alpha*2*u(1), 0;0, alpha*2*u(2)];
F_v=zeros(4,2);

Q_x=2*Q*x+fx'*V_x+F_x'*V_xx*F*1;
Q_u=2*Ru*u+fu'*V_x+F_u'*V_xx*F*1;
Q_v=-2*Rv*v+fv'*V_x+F_v'*V_xx*F*1;
Q_xx=2*Q+fx'*V_xx+V_xx*fx;
Q_ux=fu'*V_xx;
Q_uu=2*Ru;
Q_vx=fv'*V_xx;
Q_vv=-2*Rv;

dVxdt=-(Q_x-Q_ux'/Q_uu*Q_u-Q_vx'/Q_vv*Q_v);
dVxxdt=-(Q_xx-Q_ux'/Q_uu*Q_ux-Q_vx'/Q_vv*Q_vx);

dVdt=[dVxdt(:);...
    dVxxdt(:)];
end