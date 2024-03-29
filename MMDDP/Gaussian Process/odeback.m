function dVdt = odeback(t,V,u_pp,v_pp,x_pp,Q,Ru,Rv,gprMd3,gprMd4)

V_x = V(1:4);
V_xx=reshape(V(5:end), 4, 4);
u = ppval(u_pp, t);
v = ppval(v_pp, t);
x = ppval(x_pp, t);

% dmudx(1,:) = grad_gaussian(x',u',gprMd1);
% dmudx(2,:) = grad_gaussian(x',u',gprMd2);
dmudx(3,:) = grad_gaussian(x',u',gprMd3);
dmudx(4,:) = grad_gaussian(x',u',gprMd4);
fx(1,:)=[0,0,1,0];
fx(2,:)=[0,0,0,1];
fu(1,:)=[0,0];
fu(2,:)=[0,0];
for i=3:4
    fx(i,:)=dmudx(i,1:4);
    fu(i,:)=dmudx(i,5:6);
end
fv=[0,0;...
    0,0;...
    1,0;...
    0,1];

% m=1;
% l=0.5;
% g=9.81;
% FX=[                                                                                                                                                                                                                                                                                                                                                             0,                                                                                                                                                                                                                                                                                                                                                   0,                                                     1,                                                      0;...
%                                                                                                                                                                                                                                                                                                                                                              0,                                                                                                                                                                                                                                                                                                                                                   0,                                                     0,                                                      1;...
%                     - (x(3)^2*sin(x(1) - x(2))^2 - cos(x(1) - x(2))*(cos(x(1) - x(2))*x(3)^2 + x(4)^2) + (g*(2*cos(x(1)) + sin(x(1) - x(2))*sin(x(2))))/l + (u(2)*sin(x(1) - x(2)))/(l^2*m))/(cos(x(1) - x(2))^2 - 2) - (2*cos(x(1) - x(2))*sin(x(1) - x(2))*((g*(2*sin(x(1)) - cos(x(1) - x(2))*sin(x(2))))/l - sin(x(1) - x(2))*(cos(x(1) - x(2))*x(3)^2 + x(4)^2) + (u(1) - u(2)*cos(x(1) - x(2)))/(l^2*m)))/(cos(x(1) - x(2))^2 - 2)^2, (x(3)^2*sin(x(1) - x(2))^2 - cos(x(1) - x(2))*(cos(x(1) - x(2))*x(3)^2 + x(4)^2) + (g*(cos(x(1) - x(2))*cos(x(2)) + sin(x(1) - x(2))*sin(x(2))))/l + (u(2)*sin(x(1) - x(2)))/(l^2*m))/(cos(x(1) - x(2))^2 - 2) + (2*cos(x(1) - x(2))*sin(x(1) - x(2))*((g*(2*sin(x(1)) - cos(x(1) - x(2))*sin(x(2))))/l - sin(x(1) - x(2))*(cos(x(1) - x(2))*x(3)^2 + x(4)^2) + (u(1) - u(2)*cos(x(1) - x(2)))/(l^2*m)))/(cos(x(1) - x(2))^2 - 2)^2, (2*x(3)*cos(x(1) - x(2))*sin(x(1) - x(2)))/(cos(x(1) - x(2))^2 - 2),               (2*x(4)*sin(x(1) - x(2)))/(cos(x(1) - x(2))^2 - 2);...
%  - (cos(x(1) - x(2))*(2*x(3)^2 + cos(x(1) - x(2))*x(4)^2) - x(4)^2*sin(x(1) - x(2))^2 - (2*g*(cos(x(1) - x(2))*cos(x(1)) - sin(x(1) - x(2))*sin(x(1))))/l + (u(1)*sin(x(1) - x(2)))/(l^2*m))/(cos(x(1) - x(2))^2 - 2) - (2*cos(x(1) - x(2))*sin(x(1) - x(2))*(sin(x(1) - x(2))*(2*x(3)^2 + cos(x(1) - x(2))*x(4)^2) + (2*u(2) - u(1)*cos(x(1) - x(2)))/(l^2*m) + (2*g*(sin(x(2)) - cos(x(1) - x(2))*sin(x(1))))/l))/(cos(x(1) - x(2))^2 - 2)^2,      (cos(x(1) - x(2))*(2*x(3)^2 + cos(x(1) - x(2))*x(4)^2) - x(4)^2*sin(x(1) - x(2))^2 - (2*g*(cos(x(2)) - sin(x(1) - x(2))*sin(x(1))))/l + (u(1)*sin(x(1) - x(2)))/(l^2*m))/(cos(x(1) - x(2))^2 - 2) + (2*cos(x(1) - x(2))*sin(x(1) - x(2))*(sin(x(1) - x(2))*(2*x(3)^2 + cos(x(1) - x(2))*x(4)^2) + (2*u(2) - u(1)*cos(x(1) - x(2)))/(l^2*m) + (2*g*(sin(x(2)) - cos(x(1) - x(2))*sin(x(1))))/l))/(cos(x(1) - x(2))^2 - 2)^2,             -(4*x(3)*sin(x(1) - x(2)))/(cos(x(1) - x(2))^2 - 2), -(2*x(4)*cos(x(1) - x(2))*sin(x(1) - x(2)))/(cos(x(1) - x(2))^2 - 2)];
% 
% FU=[                                         0,                                         0;...
%                                          0,                                         0;...
%            -1/(l^2*m*(cos(x(1) - x(2))^2 - 2)), cos(x(1) - x(2))/(l^2*m*(cos(x(1) - x(2))^2 - 2));...
%  cos(x(1) - x(2))/(l^2*m*(cos(x(1) - x(2))^2 - 2)),           -2/(l^2*m*(cos(x(1) - x(2))^2 - 2))];
% dx=norm(FX-fx);
% du=norm(FU-fu);


Q_x=2*Q*x+fx'*V_x;
Q_u=2*Ru*u+fu'*V_x;
Q_v=-2*Rv*v+fv'*V_x;
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