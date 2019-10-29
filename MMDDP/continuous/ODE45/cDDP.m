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
N=500;
% dimenison
n = 2;  % state
m = 1;  % control
% interval


% initial
ut=linspace(0,1,N+1);
u=zeros(m,size(ut,2));
vt=linspace(0,1,N+1);
v=zeros(m,size(vt,2));

gamma=0.8;
itr=0;


fprintf('\n=========== begin Min-Max DDP ===========\n');
while(1)
    %%%% initial trajectory %%%%%
    tspan = [0 1];
    x0=[pi;0];
    [xt,x] = ode45(@(t,x) dynamics(t,x,ut,u,vt,v), tspan, x0);
    x=x';
    
%     % debug
%     plot(xt,x(1,:),xt,x(2,:),'linewidth',2);
    
    fx=zeros(2,2,size(xt,1));
    fu=zeros(2,1,size(xt,1));
    fv=zeros(2,1,size(xt,1));
    for i=1:size(xt,1)
        fx(:,:,i)=[0, 1;...
                9.81*2*cos(x(1,i)), -0.4];
        fu(:,:,i)=[0;4];
        fv(:,:,i)=[0;4];
    end
    
    tspan = [1 0];
    VN=[reshape(2*Q*x(:,size(xt,1)),[2,1]);reshape(2*Q,[4,1])];
    [Vt,V] = ode45(@(t,V) odeback(t,V,ut,u,vt,v,xt,x,fx,fu,fv,Q,Ru,Rv), tspan, VN);
    
    V_x=zeros(2,1,size(Vt,1));
    V_xx=zeros(2,2,size(Vt,1));
    for i=1:size(Vt,1)
        V_x(:,:,i)=[V(i,1);V(i,2)];
        V_xx(:,:,i)=[V(i,3),V(i,4);V(i,5),V(i,6)];
    end
    
%     % debug
%     plot(Vt,reshape(V_x(1,1,:),[1,size(V_x,3)]),Vt,reshape(V_x(2,1,:),[1,size(V_x,3)]),'linewidth',2);
    
    tspan = [0 1];
    dx0=[0;0];
    [dxt,dx] = ode45(@(t,dx) odeforw(t,dx,ut,u,vt,v,xt,fx,fu,fv,Ru,Rv,Vt,V_x,V_xx), tspan, dx0);
    dx=dx';
    
%     % debug
%     plot(dxt,dx(1,:),dxt,dx(2,:),'linewidth',2);
    
    for i=1:size(ut,2)
        for j=1:2
            for k=1:2
                Vxx(j,k)=interp1(Vt,reshape(V_xx(j,k,:),[1,size(V_xx,3)]),ut(i));
            end
        end
        for j=1:2
            Fu(j,1)=interp1(xt,reshape(fu(j,1,:),[1,size(fu,3)]),ut(i));
            Fv(j,1)=interp1(xt,reshape(fv(j,1,:),[1,size(fv,3)]),ut(i));
            Vx(j,1)=interp1(Vt,reshape(V_x(j,1,:),[1,size(V_x,3)]),ut(i));
            dX(j,1)=interp1(dxt,reshape(dx(j,:),[1,size(dxt,1)]),ut(i));
        end
        
        lu=-(2*u(1,i)*Ru+Fu'*Vx)/(2*Ru);
        Ku=-(Fu'*Vxx)/(2*Ru);
        lv=-(2*v(1,i)*Rv+Fv'*Vx)/(2*Rv);
        Kv=-(Fv'*Vxx)/(2*Rv);
        
        du(:,i)=lu+Ku*dX;
        dv(:,i)=lv+Kv*dX;
    end
    
%     % debug
%     plot(ut,du,vt,dv,'linewidth',2);
    u=u+du*gamma;
    v=v+dv*gamma;

    itr=itr+1;

%     if max(abs(du))+ max(abs(dv))< 1e-5
    if itr>=20
        break;
    end

end


fprintf(['\n'...
    'iterative   times:   %d\n'],...
    itr);
fprintf('optimal trajectory:   %.4f\n', x(:,size(xt,1)));
% fprintf('optimal u:   %.4f\n', u);
fprintf(['\n'...
    '=========== end Min-Max DDP ===========\n']);

% %% Plot
% control sequence
figure(1);
plot(ut,u,vt,v,'linewidth',2);
title('Control sequence');
xlabel('Time in sec');
ylabel('u');

% state trajectory
figure(3)
plot(xt,x(1,:),xt,x(2,:),'linewidth',2);
title('state trajectory');

end