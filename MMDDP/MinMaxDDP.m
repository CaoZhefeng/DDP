function MinMaxDDP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical
% high-dimension
% with terminal constraint
% deterministic
% continuous case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J= u_syms^2-v_syms^2;
% Phi=(x_syms^2);
% goal state: x_expected = [0,0]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cost parameter
Q=diag([1 0.1]);
Ru=0.01;
Rv=1;
% time horizon N
N=200;
% dimenison
n = 2;  % state
m = 1;  % control
% interval
dt=0.001;


% initial
u0=zeros(m,N);
u=u0;
v0=zeros(m,N);
v=v0;
x=zeros(n,N);x(:,1)=[pi,0]';

gamma=0.8;
itr=0;

% initial trajectory
[x,dx,fx,fu,fv]=dynamics(x,u,v,0,dt);

fprintf('\n=========== begin Min-Max DDP ===========\n');
while 1
    V_x(:,N)=2*Q*x(:,N);
    V_xx(:,:,N)=2*Q;
    
    for i=N:-1:2
        c_x=0;
        c_u=2*u(:,i)*Ru;
        c_v=-2*v(:,i)*Rv;
        c_xx=0;
        c_uu=2*eye(m)*Ru;
        c_vv=-2*eye(m)*Rv;
        c_ux=0;
        c_vx=0;
        c_uv=0;
        
        Q_x(:,i)=c_x+fx(:,:,i)'*V_x(:,i);
        Q_u(:,i)=c_u+fu(:,:,i)'*V_x(:,i);
        Q_v(:,i)=c_v+fv(:,:,i)'*V_x(:,i);
        Q_xx(:,:,i)=c_xx+fx(:,:,i)'*V_xx(:,:,i)+V_xx(:,:,i)*fx(:,:,i);
        Q_uu(:,:,i)=c_uu;
        Q_vv(:,:,i)=c_vv;
        Q_ux(:,:,i)=c_ux+fu(:,:,i)'*V_xx(:,:,i);
        Q_xu(:,:,i)=Q_ux(:,:,i)';
        Q_vx(:,:,i)=c_vx+fv(:,:,i)'*V_xx(:,:,i);
        Q_xv(:,:,i)=Q_vx(:,:,i)';
        Q_uv(:,:,i)=c_uv;
        Q_vu(:,:,i)=Q_uv(:,:,i)';
        
        lu(:,i)=-(Q_u(:,i)-Q_uv(:,:,i)/Q_vv(:,:,i)*Q_v(:,i))/(Q_uu(:,:,i)-Q_uv(:,:,i)/Q_vv(:,:,i)*Q_vu(:,:,i));
        lv(:,i)=-(Q_v(:,i)-Q_vu(:,:,i)/Q_uu(:,:,i)*Q_u(:,i))/(Q_vv(:,:,i)-Q_vu(:,:,i)/Q_uu(:,:,i)*Q_uv(:,:,i));
        Ku(:,:,i)=-(Q_ux(:,:,i)-Q_uv(:,:,i)/Q_vv(:,:,i)*Q_vx(:,:,i))/(Q_uu(:,:,i)-Q_uv(:,:,i)/Q_vv(:,:,i)*Q_vu(:,:,i));
        Kv(:,:,i)=-(Q_vx(:,:,i)-Q_vu(:,:,i)/Q_uu(:,:,i)*Q_ux(:,:,i))/(Q_vv(:,:,i)-Q_vu(:,:,i)/Q_uu(:,:,i)*Q_uv(:,:,i));
        
        V_x(:,i-1)=V_x(:,i)+dt*(Q_x(:,i)+Ku(:,:,i)'*Q_u(:,i)+Kv(:,:,i)'*Q_v(:,i)+Q_ux(:,:,i)'*lu(:,i)+Q_vx(:,:,i)'*lv(:,i)+...
            Ku(:,:,i)'*Q_uu(:,:,i)*lu(:,i)+Ku(:,:,i)'*Q_uv(:,:,i)*lv(:,i)+Kv(:,:,i)'*Q_vu(:,:,i)*lu(:,i)+Kv(:,:,i)'*Q_vv(:,:,i)*lv(:,i));
        V_xx(:,:,i-1)=V_xx(:,:,i)+dt*(Q_xx(:,:,i)+Ku(:,:,i)'*Q_ux(:,:,i)+Q_ux(:,:,i)'*Ku(:,:,i)+Kv(:,:,i)'*Q_vx(:,i)+Q_vx(:,i)'*Kv(:,:,i)+...
            Kv(:,:,i)'*Q_vu(:,:,i)*Ku(:,:,i)+Ku(:,:,i)'*Q_uv(:,:,i)*Kv(:,:,i)+Ku(:,:,i)'*Q_uu(:,:,i)*Ku(:,:,i)+Kv(:,:,i)'*Q_vv(:,:,i)*Ku(:,:,i));
    end
        c_x=0;
        c_u=2*u(:,1)*Ru;
        c_v=-2*v(:,1)*Rv;
        c_xx=0;
        c_uu=2*eye(m)*Ru;
        c_vv=-2*eye(m)*Rv;
        c_ux=0;
        c_vx=0;
        c_uv=0;
        
        Q_x(:,1)=c_x+fx(:,:,1)'*V_x(:,1);
        Q_u(:,1)=c_u+fu(:,:,1)'*V_x(:,1);
        Q_v(:,1)=c_v+fv(:,:,1)'*V_x(:,1);
        Q_xx(:,:,1)=c_xx+fx(:,:,1)'*V_xx(:,:,1)+V_xx(:,:,1)*fx(:,:,1);
        Q_uu(:,:,1)=c_uu;
        Q_vv(:,:,1)=c_vv;
        Q_ux(:,:,1)=c_ux+fu(:,:,1)'*V_xx(:,:,1);
        Q_xu(:,:,1)=Q_ux(:,:,1)';
        Q_vx(:,:,1)=c_vx+fv(:,:,1)'*V_xx(:,:,1);
        Q_xv(:,:,1)=Q_vx(:,:,1)';
        Q_uv(:,:,1)=c_uv;
        Q_vu(:,:,1)=Q_uv(:,:,1)';
        
        lu(:,1)=-(Q_u(:,1)-Q_uv(:,:,1)/Q_vv(:,:,1)*Q_v(:,1))/(Q_uu(:,:,1)-Q_uv(:,:,1)/Q_vv(:,:,1)*Q_vu(:,:,1));
        lv(:,1)=-(Q_v(:,1)-Q_vu(:,:,1)/Q_uu(:,:,1)*Q_u(:,1))/(Q_vv(:,:,1)-Q_vu(:,:,1)/Q_uu(:,:,1)*Q_uv(:,:,1));
        Ku(:,:,1)=-(Q_ux(:,:,1)-Q_uv(:,:,1)/Q_vv(:,:,1)*Q_vx(:,:,1))/(Q_uu(:,:,1)-Q_uv(:,:,1)/Q_vv(:,:,1)*Q_vu(:,:,1));
        Kv(:,:,1)=-(Q_vx(:,:,1)-Q_vu(:,:,1)/Q_uu(:,:,1)*Q_ux(:,:,1))/(Q_vv(:,:,1)-Q_vu(:,:,1)/Q_uu(:,:,1)*Q_uv(:,:,1));
    
    for i=1:N-1
        du(:,i)=lu(:,i)+ Ku(:,:,i)*dx(:,i);
        dv(:,i)=lv(:,i)+ Kv(:,:,i)*dx(:,i);
        u(:,i)=gamma*du(:,i)+u(:,i);
        v(:,i)=gamma*dv(:,i)+v(:,i);
        [x(:,i+1),dx(:,i+1),fx(:,:,i+1),fu(:,:,i+1),fv(:,:,i+1)]=dynamics(x,u,v,i,dt);
    end
    
    itr=itr+1;
    cost(itr)=sum(u.*u*Ru-v.*v*Rv)*dt+x(:,N)'*Q*x(:,N);
    
    max(abs(du))+max(abs(dv))
    if max(abs(du))+max(abs(dv))< 1e-4
%     if itr>=20
        break;
    end

end


fprintf(['\n'...
    'iterative   times:   %d\n'],...
    itr);
fprintf('optimal trajectory:   %.4f\n', x(:,N));
% fprintf('optimal u:   %.4f\n', u);
fprintf(['\n'...
    '=========== end BSDDP ===========\n']);

%% Plot

% control sequence
figure(1);
hold on;
% for i=1:N
%   mean_true(i)=sum(x_true(i,:))/SAM;
%   cov_true(i)=sum((x_true(i,:)-mean_true(i)*ones(1,SAM)).^2)/SAM;
% end
plot(0:dt:dt*(N-2),v(1,:),'b','linewidth',2);
plot(0:dt:dt*(N-2),u(1,:),'g','linewidth',2);

title('Control sequence');
xlabel('Time in sec');
hold off;

% Cost
figure(2);
plot(1:itr,cost,'r','linewidth',2);
title('Control sequence');
xlabel('Iteration');
end