function ddxdt = odeforw(t,dx,u_pp,v_pp,x_pp,Ru,Rv,V_pp,alpha)

u = ppval(u_pp, t);
v = ppval(v_pp, t);
x = ppval(x_pp, t);
V = ppval(V_pp, t);

V_x = V(1:2);
V_xx=reshape(V(3:end), 2, 2);

fx=[0, 1;...
    9.81*2*cos(x(1)), -0.4];
fu=[0;4];
fv=[0;4];
F=[0;alpha*u];
F_u=[0; alpha];
F_v=[0; 0];

% stochastic case
Q_u=2*Ru*u+fu'*V_x+F_u'*V_xx*F*1;
Q_v=-2*Rv*v+fv'*V_x+F_v'*V_xx*F*1;
Q_ux=fu'*V_xx;
Q_uu=2*Ru;
Q_vx=fv'*V_xx;
Q_vv=-2*Rv;

lu=-Q_u/Q_uu;
Ku=-Q_ux/Q_uu;
lv=-Q_v/Q_vv;
Kv=-Q_vx/Q_vv;

ddxdt=fu*lu+fv*lv+(fx+fu*Ku+fv*Kv)*dx;

end