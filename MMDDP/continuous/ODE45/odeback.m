function dVdt = odeback(t,V,ut,u,vt,v,xt,X,Fx,Fu,Fv,Q,Ru,Rv)
dVdt=zeros(6,1);

Vx=[V(1);V(2)];
Vxx=[V(3),V(5);V(4),V(6)];
u = interp1(ut,u,t, 'spline'); 
v = interp1(vt,v,t, 'spline'); 
for i=1:2
    for j=1:2
        fx(i,j)=interp1(xt,reshape(Fx(i,j,:),[1,size(Fx,3)]),t);
    end
end
for i=1:2
    x(i,1)=interp1(xt,reshape(X(i,:),[1,size(X,2)]),t);
    fu(i,1)=interp1(xt,reshape(Fu(i,1,:),[1,size(Fu,3)]),t);
    fv(i,1)=interp1(xt,reshape(Fv(i,1,:),[1,size(Fv,3)]),t);
end

dVxdt=-(2*Q*x+fx'*Vx-Vxx'*fu*(2*Ru*u+fu'*Vx)/(2*Ru)-Vxx'*fv*(2*Rv*v+fv'*Vx)/(2*Rv));
dVxxdt=-(2*Q+fx'*Vxx+Vxx*fx-Vxx'*(fu*fu')*Vxx/(2*Ru)-Vxx'*(fv*fv')*Vxx/(2*Rv));

dVdt(1)=dVxdt(1);
dVdt(2)=dVxdt(2);
dVdt(3)=dVxxdt(1,1);
dVdt(4)=dVxxdt(2,1);
dVdt(5)=dVxxdt(1,2);
dVdt(6)=dVxxdt(2,2);
end