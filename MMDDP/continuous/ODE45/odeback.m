function dVdt = odeback(t,V,u_pp,v_pp,x_pp,Q,Ru,Rv)

Vx = V(1:2);
Vxx=reshape(V(3:end), 2, 2);
u = ppval(u_pp, t);
v = ppval(v_pp, t);
x = ppval(x_pp, t);

fx=[0, 1;...
    9.81*2*cos(x(1)), -0.4];
fu=[0;4];
fv=[0;4];

dVxdt=-(2*Q*x+fx'*Vx-Vxx'*fu*(2*Ru*u+fu'*Vx)/(2*Ru)-Vxx'*fv*(2*Rv*v+fv'*Vx)/(2*Rv));
dVxxdt=-(2*Q+fx'*Vxx+Vxx*fx-Vxx'*(fu*fu')*Vxx/(2*Ru)-Vxx'*(fv*fv')*Vxx/(2*Rv));

dVdt=[dVxdt(:);...
    dVxxdt(:)];
end