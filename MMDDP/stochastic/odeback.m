function dVdt = odeback(t,V,u_pp,v_pp,x_pp,Q,Ru,Rv,alpha)

V_x = V(1:2);
V_xx=reshape(V(3:end), 2, 2);
u = ppval(u_pp, t);
v = ppval(v_pp, t);
x = ppval(x_pp, t);

fx=[0, 1;...
    9.81*2*cos(x(1)), -0.4];
fu=[0;4];
fv=[0;4];
F=[0;alpha*u];
F_x=[0,0;0,0];
F_u=[0; alpha];
F_v=[0; 0];

% stochastic case
% Q_x=2*Q*x+fx'*V_x+V_xx*f+F_x'*V_xx*F*1;
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