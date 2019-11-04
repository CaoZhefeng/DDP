function dDDP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical
% high-dimension
% with terminal constraint
% deterministic
% discrete case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J= u_syms^2-v_syms^2;
% Phi=(x_syms^2);
% goal state: x_expected = [0,0]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cost parameter
Q=diag([1 0.1]);
Ru=0.01;
Rv=0.8;
% time horizon N
N=100;
% dimenison
n = 2;  % state
m = 1;  % control
% interval
dt=0.01;


% initial
u0=zeros(m,N-1);
u=u0;
v0=zeros(m,N-1);
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
%     V_x(:,N)=zeros(2,1);
%     V_xx(:,:,N)=zeros(2,2);
    
    for i=N-1:-1:1
        c_x=2*Q*x(:,i)*dt;
        c_u=2*u(:,i)*Ru*dt;
        c_v=-2*v(:,i)*Rv*dt;
        c_xx=2*Q*dt;
        c_uu=2*eye(m)*Ru*dt;
        c_vv=-2*eye(m)*Rv*dt;
        c_ux=0;
        c_vx=0;
        c_uv=0;
        
        Q_x(:,i)=c_x+fx(:,:,i)'*V_x(:,i+1);
        Q_u(:,i)=c_u+fu(:,:,i)'*V_x(:,i+1);
        Q_v(:,i)=c_v+fv(:,:,i)'*V_x(:,i+1);
        Q_xx(:,:,i)=c_xx+fx(:,:,i)'*V_xx(:,:,i+1)*fx(:,:,i);
        Q_uu(:,:,i)=c_uu+fu(:,:,i)'*V_xx(:,:,i+1)*fu(:,:,i);
        Q_vv(:,:,i)=c_vv+fv(:,:,i)'*V_xx(:,:,i+1)*fv(:,:,i);
        Q_ux(:,:,i)=c_ux+fu(:,:,i)'*V_xx(:,:,i+1)*fx(:,:,i);
        Q_xu(:,:,i)=Q_ux(:,:,i)';
        Q_vx(:,:,i)=c_vx+fv(:,:,i)'*V_xx(:,:,i+1)*fx(:,:,i);
        Q_xv(:,:,i)=Q_vx(:,:,i)';
        Q_uv(:,:,i)=c_uv+fu(:,:,i)'*V_xx(:,:,i+1)*fv(:,:,i);
        Q_vu(:,:,i)=Q_uv(:,:,i)';
        
        lu(:,i)=-(Q_u(:,i)-Q_uv(:,:,i)/Q_vv(:,:,i)*Q_v(:,i))/(Q_uu(:,:,i)-Q_uv(:,:,i)/Q_vv(:,:,i)*Q_vu(:,:,i));
        lv(:,i)=-(Q_v(:,i)-Q_vu(:,:,i)/Q_uu(:,:,i)*Q_u(:,i))/(Q_vv(:,:,i)-Q_vu(:,:,i)/Q_uu(:,:,i)*Q_uv(:,:,i));
        Ku(:,:,i)=-(Q_ux(:,:,i)-Q_uv(:,:,i)/Q_vv(:,:,i)*Q_vx(:,:,i))/(Q_uu(:,:,i)-Q_uv(:,:,i)/Q_vv(:,:,i)*Q_vu(:,:,i));
        Kv(:,:,i)=-(Q_vx(:,:,i)-Q_vu(:,:,i)/Q_uu(:,:,i)*Q_ux(:,:,i))/(Q_vv(:,:,i)-Q_vu(:,:,i)/Q_uu(:,:,i)*Q_uv(:,:,i));
        
        V_x(:,i)=Q_x(:,i)+Ku(:,:,i)'*Q_u(:,i)+Kv(:,:,i)'*Q_v(:,i)+Q_ux(:,:,i)'*lu(:,i)+Q_vx(:,:,i)'*lv(:,i)+...
            Ku(:,:,i)'*Q_uu(:,:,i)*lu(:,i)+Ku(:,:,i)'*Q_uv(:,:,i)*lv(:,i)+Kv(:,:,i)'*Q_vu(:,:,i)*lu(:,i)+Kv(:,:,i)'*Q_vv(:,:,i)*lv(:,i);
        V_xx(:,:,i)=Q_xx(:,:,i)+Ku(:,:,i)'*Q_ux(:,:,i)+Q_ux(:,:,i)'*Ku(:,:,i)+Kv(:,:,i)'*Q_vx(:,i)+Q_vx(:,i)'*Kv(:,:,i)+...
            Kv(:,:,i)'*Q_vu(:,:,i)*Ku(:,:,i)+Ku(:,:,i)'*Q_uv(:,:,i)*Kv(:,:,i)+Ku(:,:,i)'*Q_uu(:,:,i)*Ku(:,:,i)+Kv(:,:,i)'*Q_vv(:,:,i)*Ku(:,:,i);
    end
    
    for i=1:N-1
        du(:,i)=lu(:,i)+ Ku(:,:,i)*dx(:,i);
        dv(:,i)=lv(:,i)+ Kv(:,:,i)*dx(:,i);
        u(:,i)=gamma*du(:,i)+u(:,i);
        v(:,i)=gamma*dv(:,i)+v(:,i);
        [x(:,i+1),dx(:,i+1),fx(:,:,i),fu(:,:,i),fv(:,:,i)]=dynamics(x,u,v,i,dt);
    end
    
    itr=itr+1;
    cost(itr)=sum(u.*u*Ru-v.*v*Rv)*dt+x(:,N)'*Q*x(:,N);
    for i=1:N-1
        cost(itr)=cost(itr)+x(:,i)'*Q*x(:,i)*dt;
    end
%     max(abs(du))+max(abs(dv))
    if max(abs(du))+max(abs(dv))< 1e-2
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
subplot(1,2,1)
% for i=1:N
%   mean_true(i)=sum(x_true(i,:))/SAM;
%   cov_true(i)=sum((x_true(i,:)-mean_true(i)*ones(1,SAM)).^2)/SAM;
% end
plot(0:dt:dt*(N-2),u(1,:),'g','linewidth',2);
title('Control sequence of u');
xlabel('Time in sec');
ylabel('u');
subplot(1,2,2)
plot(0:dt:dt*(N-2),v(1,:),'b','linewidth',2);
title('Control sequence of v');
xlabel('Time in sec');
ylabel('v');
% Cost u
figure(2);
plot(1:itr,cost,'r','linewidth',2);
title('Cost');
xlabel('Iteration');

% state trajectory
figure(3)
plot(0:dt:dt*(N-1),x(1,:),0:dt:dt*(N-1),x(2,:),'linewidth',2);
title('state trajectory');
end