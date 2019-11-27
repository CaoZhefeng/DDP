function GDDP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical
% double pendulum
% deterministic\continuous case (ODE method)
% Gaussian Process regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model-free
% Phi=(x_syms^2);
% initial state: x0 = [pi,0]';
% goal state: x_expected = [0,0]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter
% cost parameter
Q=diag([1 1 0.1 0.1]);
Ru=diag([0.01 0.01]);
Rv=diag([1 1]);
N=100; % time horizon N
n = 4;  % state dimenison
m = 2;  % control dimenison
dt=0.01; % interval
gamma=0.3;

lmd=0;

% initialize control sequence
time=0:dt:N*dt;
time_u = time(1:end-1);

u=zeros(m,N);
u_pp=spline(time_u, u);
v=zeros(2,N);
v_pp=spline(time_u, v);


%% Gaussian Process regression
[gprMd1,gprMd2,gprMd3,gprMd4]=gaussian_process;


%% Min-Max DDP
fprintf('\n=========== begin Min-Max DDP ===========\n');

itr=1;
while(1)
    %% Initial Trajectory 
    x0=[pi;pi;0;0];
    [~,x] = ode45(@(t,x) dynamics(t,x,u_pp,v_pp,gprMd1,gprMd2,gprMd3,gprMd4), time, x0);
    
    x_pp=spline(time, x');
    
    %% cost function
    terminal_cost(itr) = x(end,:)* Q * x(end,:)';
    [~, RunningCost] = ode45(@(t,J) running_cost_minmax(t, J, x_pp, u_pp, v_pp, Q, Ru, Rv), time, 0);
    cost_traj(itr) = RunningCost(end) + terminal_cost(itr);
        
    %% Backward Pass
    time_back=N*dt:-dt:0;
%     VN=[0;0;0;0;0;0];
    VN=2*[Q*x(end,:)';reshape(Q, 16, 1)];
    [~,V] = ode45(@(t,V) odeback(t,V,u_pp,v_pp,x_pp,Q,Ru,Rv,gprMd1,gprMd2,gprMd3,gprMd4), time_back, VN);
    V_pp = spline(time, V(end:-1:1, :)');
    
    V_x = V(end:-1:1, 1:4)';
    V_xx = zeros(4, 4, N+1);
    for i = 1:N+1
        V_xx(:,:,i) = reshape(V(end+1-i, 5:end), 4, 4);
    end
    
    %% Trajectory Update
    dx0=[0;0;0;0];
    [~,dx] = ode45(@(t,dx) odeforw(t,dx,u_pp,v_pp,x_pp,Ru,Rv,V_pp,gprMd1,gprMd2,gprMd3,gprMd4), time, dx0);
    
    %% Control Update
    [u, v,du_max,dv_max,lu,Ku,lv,Kv] = du_update(dx, x, u, v, V_x, V_xx, Ru,Rv,gamma,N,gprMd1,gprMd2,gprMd3,gprMd4);
    u_pp=spline(time_u, u);
    v_pp=spline(time_u, v);

    %% Convergence Check
%     if abs(x(end,1))+abs(x(end,2))< 1e-3
    if du_max+dv_max< 1e-3
%     if itr>=20
        break;
    end
    
    itr=itr+1;
end

figure(5);
plot(1:itr, cost_traj,'linewidth',2)
grid on

x=x';
fprintf(['\n'...
    'iterative   times:   %d\n'],...
    itr);
fprintf('optimal trajectory:   %.4f\n', x(:,N+1));
fprintf(['\n'...
    '=========== end Min-Max DDP ===========\n']);

% %% Stochastic test sample
% figure(4);
% hold on;
% m=1;
% l=0.5;
% g=9.81;
% alpha=0.4;
% 
% for k=1:1
%     x_true(:,1)=[pi;pi;0;0];
%     cost_sam(k)=0;
%     for i=1:N
%         dx=x_true(:,i)-x(:,i);
%         du(:,i)=lu(:,:,i)+Ku(:,:,i)*dx;
%         dv(:,i)=lv(:,:,i)+Kv(:,:,i)*dx;
% 
%         u_true(:,i)=u(:,i)+du(:,i)*gamma;
%         v_true(:,i)=v(:,i)+dv(:,i)*gamma;
%         x_true(:,i+1)=x_true(:,i)+dt*[x_true(3,i);x_true(4,i);...
%         (((u_true(1,i)+v_true(1,i))-(u_true(2,i)+v_true(2,i))*cos(x_true(1,i)-x_true(2,i)))/(m*l^2)-sin(x_true(1,i)-x_true(2,i))*(x_true(3,i)^2*cos(x_true(1,i)-x_true(2,i))+x_true(4,i)^2)+g/l*(2*sin(x_true(1,i))-sin(x_true(2,i))*cos(x_true(1,i)-x_true(2,i))))/(2-cos(x_true(1,i)-x_true(2,i))^2);...
%         ((2*(u_true(2,i)+v_true(2,i))-(u_true(1,i)+v_true(1,i))*cos(x_true(1,i)-x_true(2,i)))/(m*l^2)+sin(x_true(1,i)-x_true(2,i))*(2*x_true(3,i)^2+x_true(4,i)^2*cos(x_true(1,i)-x_true(2,i)))+2*g/l*(sin(x_true(2,i))-sin(x_true(1,i))*cos(x_true(1,i)-x_true(2,i))))/(2-cos(x_true(1,i)-x_true(2,i))^2)]+...
%         [0; 0; alpha*u_true(1,i)*0.1*randn; alpha*u_true(2,i)*0.1*randn];
%         cost_sam(k)=cost_sam(k)+dt*(x_true(:,i)'* Q * x_true(:,i) + u_true(:,i)' * Ru * u_true(:,i) - v_true(:,i)' * Rv * v_true(:,i));
%     end
%     cost_sam(k)=cost_sam(k)+x_true(:,N+1)'* Q * x_true(:,N+1);
%     plot(time,x_true(1,:),time,x_true(2,:),time,x_true(3,:),time,x_true(4,:),'linewidth',2);
% end
% cost=mean(cost_sam)

%% Plot
% control sequence
figure(1);
plot(time_u,u(1,:),time_u,u(2,:),'linewidth',2);
legend('Tau1','Tau2');

figure(2);
plot(time_u,v(1,:),time_u,v(2,:),'linewidth',2);
legend('v1','v2');
% title('Control sequence');
% xlabel('Time in sec');
% ylabel('u\\v');

% state trajectory
figure(3)
% plot(time,x(1,:),time,x(2,:),'linewidth',2);
plot(time,x(1,:),time,x(2,:),time,x(3,:),time,x(4,:),'linewidth',2);
legend('theta1','theta2','omega1','omega2');
title('state trajectory');

end