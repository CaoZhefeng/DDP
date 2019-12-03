function [unew, vnew,du_max,dv_max,lu,Ku,lv,Kv] = du_update(dx, X, u, v, V_x, V_xx, Ru,Rv,gamma,N)
unew=u*0;
vnew=v*0;
m=1;
l=0.5;
g=9.81;

for i=1:N
    x=reshape(X(i,:),4,1);
    fu=[                                         0,                                         0;...
                                         0,                                         0;...
           -1/(l^2*m*(cos(x(1) - x(2))^2 - 2)), cos(x(1) - x(2))/(l^2*m*(cos(x(1) - x(2))^2 - 2));...
 cos(x(1) - x(2))/(l^2*m*(cos(x(1) - x(2))^2 - 2)),           -2/(l^2*m*(cos(x(1) - x(2))^2 - 2))];
    fv=[0,0;...
        0,0;...
        1,0;...
         0,1];
    
    Q_u=2*Ru*u(:,i)+fu'*V_x(:,i);
    Q_v=-2*Rv*v(:,i)+fv'*V_x(:,i);
    Q_ux=fu'*V_xx(:,:,i);
    Q_uu=2*Ru;
    Q_vx=fv'*V_xx(:,:,i);
    Q_vv=-2*Rv;

    lu(:,:,i)=-Q_uu\Q_u;
    Ku(:,:,i)=-Q_uu\Q_ux;
    lv(:,:,i)=-Q_vv\Q_v;
    Kv(:,:,i)=-Q_vv\Q_vx;

    du(:,i)=lu(:,:,i)+Ku(:,:,i)*dx(i,:)';
    dv(:,i)=lv(:,:,i)+Kv(:,:,i)*dx(i,:)';
    
    unew(:,i)=u(:,i)+du(:,i)*gamma;
    vnew(:,i)=v(:,i)+dv(:,i)*gamma;
end
du_max=max(abs(du(:)));
dv_max=max(abs(dv(:)));

end