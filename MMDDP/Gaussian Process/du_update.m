function [unew, vnew,du_max,dv_max,lu,Ku,lv,Kv] = du_update(dx, X, u, v, V_x, V_xx, Ru,Rv,gamma,N,gprMd3,gprMd4)
unew=u*0;
vnew=v*0;

for i=1:N
    x=reshape(X(i,:),4,1);
%     dmudx(1,:) = grad_gaussian(x',u(:,i)',gprMd1);
%     dmudx(2,:) = grad_gaussian(x',u(:,i)',gprMd2);
    dmudx(3,:) = grad_gaussian(x',u(:,i)',gprMd3);
    dmudx(4,:) = grad_gaussian(x',u(:,i)',gprMd4);
    fu(1,:)=[0,0];
    fu(2,:)=[0,0];
    for j=3:4
        fu(j,:)=dmudx(j,5:6);
    end
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