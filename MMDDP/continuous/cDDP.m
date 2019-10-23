function cDDP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical
% high-dimension
% with terminal constraint
% deterministic
% continuous case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J= u^2-v^2+x^2;
% Phi=(x_syms^2);
% goal state: x_expected = [0,0]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cost parameter
Q=diag([1 0.1]);
Ru=0.01;

% time horizon N
N=100;
% dimenison
n = 2;  % state
m = 1;  % control
% interval
dt=0.01;

% initial
u=zeros(m,N);
x=zeros(n,N);x(:,1)=[pi,0]';
dx=zeros(n,N);dx(:,1)=[0;0];
gamma=0.8;
itr=0;


fprintf('\n=========== begin Min-Max DDP ===========\n');
while(1)
    %%%% initial trajectory %%%%%
    [x, fx, fu]=dynamics(x,u,dt);
    % tspan = [0 1];
    % x0=[pi;0];
    % [xt,x] = ode45(@(t,x) dynamics(t,x,ut,u), tspan, x0);
    % for i=1:size(xt)
    %     fx(:,:,i)=[0, 1;...
    %             9.81*2*cos(x(i,1)), -0.4];
    % end

    % calculate value function at final time
    V_x(:,N)=2*Q*x(:,N);
    V_xx(:,:,N)=2*Q;

    for i=N-1:-1:1
        % Q_x\ Q_xx change here
        V_x(:,i)=V_x(:,i+1)+dt*(2*Q*x(:,i+1)+fx(:,:,i+1)'*V_x(:,i+1)-V_xx(:,:,i+1)'*fu(:,:,i+1)*(2*Ru*u(:,i)+fu(:,:,i+1)'*V_x(:,i+1))/(2*Ru));
        V_xx(:,:,i)=V_xx(:,:,i+1)+dt*(2*Q+V_xx(:,:,i+1)*fx(:,:,i+1)+fx(:,:,i+1)'*V_xx(:,:,i+1)-V_xx(:,:,i+1)'*fu(:,:,i+1)*fu(:,:,i+1)'*V_xx(:,:,i+1)/(2*Ru));

        Q_u(:,i)=2*u(:,i)*Ru+fu(:,:,i)'*V_x(:,i);
        Q_uu(:,:,i)=2*eye(m)*Ru;
        Q_ux(:,:,i)=fu(:,:,i)'*V_xx(:,:,i+1);
        lu(:,i)=-Q_u(:,i)/Q_uu(:,:,i);
        Ku(:,:,i)=-Q_ux(:,:,i)/Q_uu(:,:,i); 
    end

    for i=1:N-1
        dx(:,i+1)=dx(:,i)+dt*(fu(:,:,i)*lu(:,i)+(fx(:,:,i)+fu(:,:,i)*Ku(:,:,i))*dx(:,i));
        du(:,i)=lu(:,i)+Ku(:,:,i)*dx(:,i);
    end
    du(:,N)=0;
    u=u+du*gamma;

    itr=itr+1;
%     cost(itr)=sum(u.*u*Ru)*dt+x(:,N)'*Q*x(:,N);

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

% %% Plot
% 
% % control sequence
% figure(1);
% hold on;
% % for i=1:N
% %   mean_true(i)=sum(x_true(i,:))/SAM;
% %   cov_true(i)=sum((x_true(i,:)-mean_true(i)*ones(1,SAM)).^2)/SAM;
% % end
% plot(0:dt:dt*(N-2),v(1,:),'b','linewidth',2);
% plot(0:dt:dt*(N-2),u(1,:),'g','linewidth',2);
% 
% title('Control sequence');
% xlabel('Time in sec');
% hold off;
% 
% % Cost
% figure(2);
% plot(1:itr,cost,'r','linewidth',2);
% title('Control sequence');
% xlabel('Iteration');
% end