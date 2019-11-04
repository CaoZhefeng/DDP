function [unew, vnew,du_max,dv_max,lu,Ku,lv,Kv] = du_update(dx, u, v, V_x, V_xx, Ru,Rv,gamma,N)
unew=u*0;
vnew=v*0;

for i=1:N
    
    fu=[0;4];
    fv=[0;4];
    lu(i)=-(2*u(:,i)*Ru+fu'*V_x(:,i))/(2*Ru);
    Ku(:,:,i)=-(fu'*V_xx(:,:,i))/(2*Ru);
    lv(i)=-(2*v(:,i)*Rv+fv'*V_x(:,i))/(2*Rv);
    Kv(:,:,i)=-(fv'*V_xx(:,:,i))/(2*Rv);

    du(:,i)=lu(i)+Ku(:,:,i)*dx(i,:)';
    dv(:,i)=lv(i)+Kv(:,:,i)*dx(i,:)';
    
    unew(:,i)=u(:,i)+du(:,i)*gamma;
    vnew(:,i)=v(:,i)+dv(:,i)*gamma;
end
du_max=max(abs(du));
dv_max=max(abs(dv));

end