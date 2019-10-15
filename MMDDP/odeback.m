function dVdt = odeback(t,y,ut,u,vt,v,fxt,fx,fut,fu,fvt,fv,Ru,Rv)
% u,v 是关于t的函数
Vx=[y(1);y(2)];
Vxx=[y(3),y(4);y(5),y(6)];
u = interp1(ut,u,t, 'spline'); 
v = interp1(vt,v,t, 'spline'); 
fx = interp1(fxt,fx,t, 'spline'); 
fu = interp1(fut,fu,t, 'spline'); 
fv = interp1(fvt,fv,t, 'spline'); 

dVxdt=fx'*Vx-Vxx'*fu*(2*Ru*u+fu'*Vx)/(2*Ru)-Vxx'*fv*(2*Rv*v+fv'*Vx)/(2*Rv);
Vxxdt=fx'*Vxx+Vxx*fx-Vxx'*(fu*fu')*Vxx/(2*Ru)-Vxx'*(fv*fv')*Vxx/(2*Rv);
dVdt(1)=dVxdt(1);
dVdt(2)=dVxdt(2);
dVdt(3)=dVxxdt(1,1);
dVdt(4)=dVxxdt(1,2);
dVdt(5)=dVxxdt(2,1);
dVdt(6)=dVxxdt(2,2);
end