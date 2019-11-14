function SDDP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical
% pendulum
% Stochastic\continuous case
% ODE method: ODE45
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J= x^2+u^2-v^2+x^2;
% Phi=(x_syms^2);
% initial state: x0 = [pi,0]';
% goal state: x_expected = [0,0]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter
% cost parameter
Q=diag([1 0.1]);
Ru=0.01;
Rv=1;
alpha=0.1; % noise
N=100; % time horizon N
n = 2;  % state dimenison
m = 1;  % control dimenison
dt=0.01; % interval
gamma=0.8;

% initialize control sequence
time=0:dt:N*dt;
time_u = time(1:end-1);

u=zeros(m,N);
u_pp=spline(time_u, u);
v=zeros(m,N);
v_pp=spline(time_u, v);

%% Min-Max DDP
fprintf('\n=========== begin Min-Max DDP ===========\n');

itr=1;
while(1)
    %% Initial Trajectory 
    x0=[pi;0];
    [~,x] = ode45(@(t,x) dynamics(t,x,u_pp,v_pp), time, x0);
    
    x_pp=spline(time, x');
    
    %% Backward Pass
    time_back=N*dt:-dt:0;
    VN=[reshape(2*Q*x(end,:)',[2,1]);reshape(2*Q,[4,1])];
%     VN=[0;0;0;0;0;0];
    [~,V] = ode45(@(t,V) odeback(t,V,u_pp,v_pp,x_pp,Q,Ru,Rv,alpha), time_back, VN);
    V_pp = spline(time, V(end:-1:1, :)');
    
    V_x = V(end:-1:1, 1:2)';
    V_xx = zeros(2, 2, N+1);
    for i = 1:N+1
        V_xx(:,:,i) = reshape(V(end+1-i, 3:end), 2, 2);
    end
    
    %% Trajectory Update
    dx0=[0;0];
    [~,dx] = ode45(@(t,dx) odeforw(t,dx,u_pp,v_pp,x_pp,Ru,Rv,V_pp,alpha), time, dx0);
    
    %% Control Update
    [u, v,du_max,dv_max,lu,Ku,lv,Kv] = du_update(dx, u, v, V_x, V_xx, Ru,Rv,gamma,N,alpha);
    u_pp=spline(time_u, u);
    v_pp=spline(time_u, v);

    %% Convergence Check
    if du_max+dv_max< 1e-5
%     if itr>=20
        break;
    end
    
    itr=itr+1;
end

x=x';
fprintf(['\n'...
    'iterative   times:   %d\n'],...
    itr);
fprintf('optimal trajectory:   %.4f\n', x(:,N+1));
fprintf(['\n'...
    '=========== end Min-Max DDP ===========\n']);


%% Stochastic test sample
figure(3);
hold on;
for sample=1:10
    x_true(:,1)=[pi;0];
    for i=1:N
        dx=x_true(:,i)-x(:,i);
        du(:,i)=lu(i)+Ku(:,:,i)*dx;
        dv(:,i)=lv(i)+Kv(:,:,i)*dx;

        u_true(:,i)=u(:,i)+du(:,i)*gamma;
        v_true(:,i)=v(:,i)+dv(:,i)*gamma;
        x_true(:,i+1)=x_true(:,i)+dt*[x_true(2,i); 9.81*2*sin(x_true(1,i))-0.4*x_true(2,i)+4*u_true(1,i)+4*v_true(1,i)]+[0; alpha* u_true(1,i)]*0.1*randn;
    end
    plot(x_true(1,:),x_true(2,:),'linewidth',2);
end

%% Plot
% control sequence
figure(1);
plot(time_u,u,time_u,v,'linewidth',2);
title('Control sequence');
xlabel('Time in sec');
ylabel('u\\v');

% state trajectory
figure(2)
plot(x(1,:),x(2,:),'linewidth',2);
title('state trajectory');

end