function [gprMd1,gprMd2,gprMd3,gprMd4]=gaussian_process
%%
% learn gaussian process for double pendulum
load nom_tra.mat

rng(10)
number_Observation = 101;
% X_doupend = [-pi+(2*pi)*rand(number_Observation,1),-pi+(2*pi)*rand(number_Observation,1),...
%     -2+(4)*rand(number_Observation,1),-2+(4)*rand(number_Observation,1)];
% X_doupend = [X_doupend;[pi,pi,0,0]];
X_doupend=x_norm(:,1:100);
X_doupend=X_doupend';
X_doupend = [X_doupend;[pi,pi,0,0]];
% u = [-5+10*rand(number_Observation,1),-5+10*rand(number_Observation,1)];
% u = [u;[0,0]];
u=u_norm';
u = [u;[0,0]];
number_Observation=number_Observation+1;
% v = 0.03*[-5+10*rand(number_Observation,1),-5+10*rand(number_Observation,1)];
dxdt_doupend = zeros(number_Observation,4);
dxdt_doupend = vectorized_dynamics(X_doupend,u);
X_doupend = [X_doupend,u];

sigma_noise_1=std(dxdt_doupend(:,1));
sigma_noise_2=std(dxdt_doupend(:,2));
sigma_noise_3=0.1*std(dxdt_doupend(:,3));
sigma_noise_4=0.1*std(dxdt_doupend(:,4));
sigmaM0 = 0.8*ones(6,1);


gprMd1 = fitrgp(X_doupend,dxdt_doupend(:,1),'Basis','linear','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigma_noise_1],'Sigma',sigma_noise_1,'Standardize',1);

gprMd2 = fitrgp(X_doupend,dxdt_doupend(:,2),'Basis','linear','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigma_noise_1],'Sigma',sigma_noise_2,'Standardize',1);

gprMd3 = fitrgp(X_doupend,dxdt_doupend(:,3),'Basis','constant','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigma_noise_1],'Sigma',sigma_noise_3,'Standardize',1);

gprMd4 = fitrgp(X_doupend,dxdt_doupend(:,4),'Basis','constant','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigma_noise_1],'Sigma',sigma_noise_4,'Standardize',1);

%% predict gradient

% if model employ the constant regression method than we need to consider
% the differentiation of gaussian process otherwise, take beta into
% consideration
end