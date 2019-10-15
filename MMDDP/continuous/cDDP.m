function cDDP
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

% time horizon N
N=200;
% dimenison
n = 2;  % state
m = 1;  % control
% interval
dt=0.01;


% initial
u0=zeros(m,N);
u=u0;
x=zeros(n,N);x(:,1)=[pi,0]';

gamma=0.8;
itr=0;

% initial trajectory
[x,dx,fx,fu]=dynamics(x,u,0,dt);

fprintf('\n=========== begin Min-Max DDP ===========\n');
while 1
    V_x(:,N)=2*Q*x(:,N);
    V_xx(:,:,N)=2*Q;
    
    for i=N-1:-1:1
        c_x=0;
        c_u=2*u(:,i)*Ru*dt;
        c_xx=0;
        c_uu=2*eye(m)*Ru*dt;
        c_ux=0;
        
        Q_x(:,i)=c_x+fx(:,:,i)'*V_x(:,i+1);
        Q_u(:,i)=c_u+fu(:,:,i)'*V_x(:,i+1);
        Q_xx(:,:,i)=c_xx+fx(:,:,i)'*V_xx(:,:,i+1)+V_xx(:,:,i+1)*fx(:,:,i);
        Q_uu(:,:,i)=c_uu;
        Q_ux(:,:,i)=c_ux+fu(:,:,i)'*V_xx(:,:,i+1);
        Q_xu(:,:,i)=Q_ux(:,:,i)';
        
        lu(:,i)=-Q_u(:,i)/Q_uu(:,:,i);
        Ku(:,:,i)=-Q_ux(:,:,i)/Q_uu(:,:,i);
        
        V_x(:,i)=V_x(:,i+1)+dt*(Q_x(:,i)+Ku(:,:,i)'*Q_u(:,i)+Q_ux(:,:,i)'*lu(:,i)+Ku(:,:,i)'*Q_uu(:,:,i)*lu(:,i));
        V_xx(:,:,i)=V_xx(:,:,i+1)+dt*(Q_xx(:,:,i)+Ku(:,:,i)'*Q_ux(:,:,i)+Q_ux(:,:,i)'*Ku(:,:,i)+Ku(:,:,i)'*Q_uu(:,:,i)*Ku(:,:,i));
    end
    
    for i=1:N-1
        du(:,i)=lu(:,i)+ Ku(:,:,i)*dx(:,i);
        u(:,i)=gamma*du(:,i)+u(:,i);
        [x(:,i+1),dx(:,i+1),fx(:,:,i+1),fu(:,:,i+1)]=dynamics(x,u,i,dt);
    end
    
    itr=itr+1;
    cost(itr)=sum(u.*u*Ru)*dt+x(:,N)'*Q*x(:,N)
    
    max(abs(du))
    if max(abs(du))< 1e-5
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