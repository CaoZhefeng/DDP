function ddxdt = odeforw(t,dx,ut,u,vt,v,xt,Fx,Fu,Fv,Ru,Rv,Vt,Vx,Vxx)
ddxdt=zeros(2,1);

u = interp1(ut,u,t, 'spline'); 
v = interp1(vt,v,t, 'spline'); 
for i=1:2
    for j=1:2
        fx(i,j)=interp1(xt,reshape(Fx(i,j,:),[1,size(Fx,3)]),t);
        V_xx(i,j)=interp1(Vt,reshape(Vxx(i,j,:),[1,size(Vxx,3)]),t);
    end
end
for i=1:2
    fu(i,1)=interp1(xt,reshape(Fu(i,1,:),[1,size(Fu,3)]),t);
    fv(i,1)=interp1(xt,reshape(Fv(i,1,:),[1,size(Fv,3)]),t);
    V_x(i,1)=interp1(Vt,reshape(Vx(i,1,:),[1,size(Vx,3)]),t);
end


lu=-(2*u*Ru+fu'*V_x)/(2*Ru);
Ku=-(fu'*V_xx)/(2*Ru);
lv=-(2*v*Rv+fv'*V_x)/(2*Rv);
Kv=-(fv'*V_xx)/(2*Rv);

ddxdt=fu*lu+fv*lv+(fx+fu*Ku+fv*Kv)*dx;

end