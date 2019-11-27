
%%
% learn gaussian process for double pendulum
X_doupend = [linspace(-pi,pi,100)',linspace(-pi,pi,100)',linspace(-10,10,100)',linspace(-10,10,100)'];
u = [linspace(-5,5,100)',linspace(-5,5,100)'];
% v = [linspace(-5,5,100)',linspace(-5,5,100)'];
dxdt_doupend = zeros(100,4);
dxdt_doupend = vectorized_dynamics(X_doupend,u);

sigma_noise_1=std(dxdt_doupend(:,1));
sigma_noise_2=std(dxdt_doupend(:,2));
sigma_noise_3=std(dxdt_doupend(:,3));
sigma_noise_4=std(dxdt_doupend(:,4));
sigmaM0 = 10*ones(4,1);

X_test = [linspace(-1,1,19)',linspace(-1,1,19)',linspace(-2.6,2.6,19)',linspace(-2.6,2.6,19)'];
u = [linspace(-4.9,3.9,19)',linspace(-4.9,3.9,19)'];
% v = [linspace(3.9,4.9,19)',linspace(3.9,4.9,19)'];
y_test = vectorized_dynamics(X_test,u);
gprMd1 = fitrgp(X_doupend,dxdt_doupend(:,1),'Basis','constant','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigma_noise_1],'Sigma',sigma_noise_1,'Standardize',1);

gprMd2 = fitrgp(X_doupend,dxdt_doupend(:,2),'Basis','constant','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigma_noise_1],'Sigma',sigma_noise_2,'Standardize',1);

gprMd3 = fitrgp(X_doupend,dxdt_doupend(:,3),'Basis','constant','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigma_noise_1],'Sigma',sigma_noise_3,'Standardize',1);

gprMd4 = fitrgp(X_doupend,dxdt_doupend(:,4),'Basis','constant','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigma_noise_1],'Sigma',sigma_noise_4,'Standardize',1);

[mu_dxdt_3,sigma_dxdt_3] = predict(gprMd3,X_test)