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
Rv=1;

% time horizon N
N=1000;
% dimenison
n = 2;  % state
m = 1;  % control
% interval
dt=0.001;

% initial
u=zeros(m,N);
v=zeros(m,N);
x=zeros(n,N);x(:,1)=[pi,0]';
dx=zeros(n,N);dx(:,1)=[0;0];
gamma=0.8;
itr=0;


fprintf('\n=========== begin Min-Max DDP ===========\n');
while(1)
    %%%% initial trajectory %%%%%
    [x, fx, fu, fv]=dynamics(x,u,v,dt);
    % tspan = [0 1];
    % x0=[pi;0];
    % [xt,x] = ode45(@(t,x) dynamics(t,x,ut,u), tspan, x0);
    % for i=1:size(xt)
    %     fx(:,:,i)=[0, 1;...
    %             9.81*2*cos(x(i,1)), -0.4];
    % end

    % calculate value function at final time
%     V_x(:,N)=2*Q*x(:,N);
%     V_xx(:,:,N)=2*Q;
    V_x(:,N)=2*Q*x(:,N);
    V_xx(:,:,N)=2*Q;

    for i=N-1:-1:1
        % Q_x\ Q_xx change here
        V_x(:,i)=V_x(:,i+1)+dt*(2*Q*x(:,i+1)+fx(:,:,i+1)'*V_x(:,i+1)-V_xx(:,:,i+1)'*fu(:,:,i+1)*(2*Ru*u(:,i)+fu(:,:,i+1)'*V_x(:,i+1))/(2*Ru)-V_xx(:,:,i+1)'*fv(:,:,i+1)*(2*Rv*v(:,i)+fv(:,:,i+1)'*V_x(:,i+1))/(2*Rv));
        V_xx(:,:,i)=V_xx(:,:,i+1)+dt*(2*Q+V_xx(:,:,i+1)*fx(:,:,i+1)+fx(:,:,i+1)'*V_xx(:,:,i+1)-V_xx(:,:,i+1)'*fu(:,:,i+1)*fu(:,:,i+1)'*V_xx(:,:,i+1)/(2*Ru)-V_xx(:,:,i+1)'*fv(:,:,i+1)*fv(:,:,i+1)'*V_xx(:,:,i+1)/(2*Rv));

        Q_u(:,i)=2*u(:,i)*Ru+fu(:,:,i)'*V_x(:,i);
        Q_uu(:,:,i)=2*eye(m)*Ru;
        Q_ux(:,:,i)=fu(:,:,i)'*V_xx(:,:,i+1);
        Q_v(:,i)=-2*v(:,i)*Rv+fv(:,:,i)'*V_x(:,i);
        Q_vv(:,:,i)=-2*eye(m)*Rv;
        Q_vx(:,:,i)=fv(:,:,i)'*V_xx(:,:,i+1);
        lu(:,i)=-Q_u(:,i)/Q_uu(:,:,i);
        Ku(:,:,i)=-Q_ux(:,:,i)/Q_uu(:,:,i);
        lv(:,i)=-Q_v(:,i)/Q_vv(:,:,i);
        Kv(:,:,i)=-Q_vx(:,:,i)/Q_vv(:,:,i); 
    end

    for i=1:N-1
        dx(:,i+1)=dx(:,i)+dt*(fu(:,:,i)*lu(:,i)+fv(:,:,i)*lv(:,i)+(fx(:,:,i)+fu(:,:,i)*Ku(:,:,i)+fv(:,:,i)*Kv(:,:,i))*dx(:,i));
        du(:,i)=lu(:,i)+Ku(:,:,i)*dx(:,i);
        dv(:,i)=lv(:,i)+Kv(:,:,i)*dx(:,i);
    end
    du(:,N)=0;
    dv(:,N)=0;
    u=u+du*gamma;
    v=v+dv*gamma;

    itr=itr+1;
%     cost(itr)=sum(u.*u*Ru)*dt+x(:,N)'*Q*x(:,N);

    if max(abs(du))+ max(abs(dv))< 1e-5
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
% control sequence
figure(1);
plot(0:dt:dt*(N-1),u(1,:),0:dt:dt*(N-1),v(1,:),'linewidth',2);
title('Control sequence');
xlabel('Time in sec');
ylabel('u');

% state trajectory
figure(3)
plot(0:dt:dt*(N-1),x(1,:),0:dt:dt*(N-1),x(2,:),'linewidth',2);
title('state trajectory');

end