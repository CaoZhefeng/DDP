function ddxdt = odeforw(t,dx,u_pp,v_pp,x_pp,Ru,Rv,V_pp)

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


lu=-(2*u*Ru+fu'*V_x)/(2*Ru);
Ku=-(fu'*V_xx)/(2*Ru);
lv=-(2*v*Rv+fv'*V_x)/(2*Rv);
Kv=-(fv'*V_xx)/(2*Rv);

ddxdt=fu*lu+fv*lv+(fx+fu*Ku+fv*Kv)*dx;

end