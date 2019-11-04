function [unew, vnew,du_max,dv_max,lu,Ku,lv,Kv] = du_update(dx, u, v, V_x, V_xx, Ru,Rv,gamma,N,alpha)
unew=u*0;
vnew=v*0;

for i=1:N
    
    fu=[0;4];
    fv=[0;4];
    F_u=[0; alpha];
    F_v=[0; 0];
    F=[0;alpha*u(:,i)];
    
    % stochastic case
    Q_u=2*Ru*u(:,i)+fu'*V_x(:,i)+F_u'*V_xx(:,:,i)*F*1;
    Q_v=-2*Rv*v(:,i)+fv'*V_x(:,i)+F_v'*V_xx(:,:,i)*F*1;
    Q_ux=fu'*V_xx(:,:,i);
    Q_uu=2*Ru;
    Q_vx=fv'*V_xx(:,:,i);
    Q_vv=-2*Rv;

    lu(i)=-Q_u/Q_uu;
    Ku(:,:,i)=-Q_ux/Q_uu;
    lv(i)=-Q_v/Q_vv;
    Kv(:,:,i)=-Q_vx/Q_vv;

    du(:,i)=lu(i)+Ku(:,:,i)*dx(i,:)';
    dv(:,i)=lv(i)+Kv(:,:,i)*dx(i,:)';
    
    unew(:,i)=u(:,i)+du(:,i)*gamma;
    vnew(:,i)=v(:,i)+dv(:,i)*gamma;
end
du_max=max(abs(du));
dv_max=max(abs(dv));

end