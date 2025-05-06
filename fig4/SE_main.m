clc
clear 


% para
L=15; % number of channel paths
f_c = 2.1e9;
c = 3e8;
lambda = c/f_c;               % \lambda
A_max  = lambda*4;           % max of A
A_min  = lambda*2;           % min of A
N_A    = 5;                 % number of different A

K = 5;  %user_num
Nh      = 3;     
Nv      = 3;  
N = Nh*Nv;
% N
Mh      = 3; 
Mv     =3;
M= Mh*Mv;
%Noise figure
  noisefigure=9;
  %communication bandwidth
  B=20e6;
  noisepowerdbm=-174+10*log10(B)+noisefigure;%%noise power
  noise=db2pow(noisepowerdbm)/1000;
% SNR        = 15;             % SNR dB
% SNR_linear = 10^(SNR/10);    % SNR
% noise = 1/SNR;
p_list = -5:5:30;
p_list =db2pow(p_list)/1000;
N_P=length(p_list);
epsil      = 10^-3; 
d          = lambda/2; 

Nuser  = 10;
Nter       = 50;%%蒙特卡洛次数


SE_ROMA= zeros(N_P, Nuser);
SE_MA= zeros(N_P, Nuser);
SE_RO= zeros(N_P, Nuser);
SE_FIS= zeros(N_P, Nuser);

%%用户分布范围
x_max = 10;
x_min = -10;
y_max = 10;
y_min = -10;
z_max =10;
z_min = -10;
tic; 



for n_user = 1:Nuser
location_u = zeros(K,3);
%%用户位置
location_u(:,1)=(x_max-x_min)*rand(K,1)+x_min;
location_u(:,2)=(y_max-y_min)*rand(K,1)+y_min;
location_u(:,3)=(z_max-z_min)*rand(K,1)+z_min;
cap        = zeros(N_P, Nter);     % capacity
cap_ma = zeros(N_P, Nter);   
cap_ro = zeros(N_P,Nter);
cap_fis = zeros(N_P, Nter);   
cap_sa = zeros(N_P, Nter);  

A =  (A_max-A_min)/(N_A-1)*(5-1) + A_min;%%固定天线孔径
delta =A/Mh;%%天线间距初始 

%%需要优化的参数有天线相对平面的位置，两个旋转角度
SE_base=zeros(N_P,Nter);
for n_A = 1:N_P
p_max =p_list(n_A);


% Initialization of alpha and beta for M=N=4, d=lambda/2
alpha = rand(K+1,1)*pi;
beta = rand(K+1,1)*pi;

% R_s2 = position(0, 0, 0, 0,0,Mh2,Mv2,delta,0); % shape 2*4
% R_r2 = zeros(N2,3,K);
% for u=1:K
% R_r2(:,:,u) = position(location_u(u,1), location_u(u,2), location_u(u,3), 0,0,Nh2,Nv2,delta,0);
% end

R_s = position(0, 0, 0, 0,0,Mh,Mv,delta,0); % shape 2*4
R_r = zeros(N,3,K);
for u=1:K
R_r(:,:,u) = position(location_u(u,1), location_u(u,2), location_u(u,3), 0,0,Nh,Nv,delta,0);
end
% Initialization of Z
Z_S = (2*rand(M, 3)-1)*A;
Z_R = (2*rand(N,3,K)-1)*A;
for iter = 1:Nter

%%计算路径损耗
path = zeros(K,1);
path_loss = zeros(K,L);
for u=1:K
distances = sqrt(location_u(u,1)^2+location_u(u,2)^2+location_u(u,3));
[path(u),~] =  functionChannelgain(distances);
path(u) = sqrt(db2pow(path(u)));
z = (randn(1, L) + 1i * randn(1, L));  % 随机生成复数数列，实部和虚部均为正态分布
% 计算当前数列的2范数
current_norm = norm(z, 'fro');

path_loss(u,:) = path(u)*z/current_norm; % 生成复数序列，模为1

end
theta_t = rand(L,1)*pi;                                          % \theta_t
phi_t   = rand(L,1)*pi;                                          % \phi_t
theta_r = rand(L,1)*pi;                                          % \theta_r
phi_r   = rand(L,1)*pi;                                          % \phi_r

% H
H_channel = zeros(M,N,K);
W_precoding = zeros(M,N,K);
B=zeros(M,M);
for u =1:K
    r_t = position(0, 0, 0, alpha(1), beta(1),Mh,Mv,delta,R_s);
    r_r = position(location_u(u,1),location_u(u,2),location_u(u,3), alpha(u+1),beta(u+1),Nh,Nv,delta,R_r(:,:,u));
    H_channel(:,:,u)=Channel(L, path_loss(u,:), theta_t, phi_t, theta_r, phi_r,...
 lambda, r_t, r_r, M, N);
    B=B+H_channel(:,:,u)*H_channel(:,:,u)';
end
for u=1:K
    W_precoding(:,:,u) = inv(B)*H_channel(:,:,u);
    W_precoding(:,:,u) =W_precoding(:,:,u)*(sqrt(p_max)/norm(W_precoding(:,:,u), 'fro'));
end
SE= SE_compute(H_channel, W_precoding, K, N, noise);
SE_base(n_A,iter)=SE;
cap_temp = SE;
cap_before = -1;
ROMA_SE = ROMA_rate(p_max,d,epsil,cap_temp, cap_before,iter,Nter, n_A, noise,...
    location_u,Mh,Mv,Nh,Nv,K,L,Z_S, R_s, R_r, theta_r, phi_r, ...
    theta_t,phi_t, alpha,beta,path_loss, lambda, A,Z_R);
cap(n_A, iter) = ROMA_SE ;
 MA_SE = MA_rate(p_max,d,epsil,cap_temp, cap_before,iter,Nter, n_A, noise,...
    location_u,Mh,Mv,Nh,Nv,K,L,Z_S, R_s, R_r, theta_r, phi_r, ...
    theta_t,phi_t, zeros(K+1,1),zeros(K+1,1),path_loss, lambda, A,Z_R);
cap_ma(n_A, iter)=MA_SE;
RO_SE = RO_rate(p_max,d,epsil,cap_temp, cap_before,iter,Nter, n_A, noise,...
    location_u,Mh,Mv,Nh,Nv,K,L, R_s, R_r, theta_r, phi_r, ...
    theta_t,phi_t, alpha,beta,path_loss, lambda, A);
cap_ro(n_A, iter)=RO_SE;
FIS_SE = FIS_rate(p_max,noise,...
    location_u,Mh,Mv,Nh,Nv,K,L, R_s, R_r, theta_r, phi_r, ...
    theta_t,phi_t, zeros(K+1,1),zeros(K+1,1),path_loss, lambda);
cap_fis(n_A, iter)=FIS_SE;
b=toc;
fprintf('ROMA_SE=%f, MA_SE=%f, FIS_SE=%f, RO_SE=%f, time=%d\n',ROMA_SE,MA_SE,FIS_SE,RO_SE,b);
end

end
SE_ROMA(:,n_user)= abs(mean(cap,2));
SE_MA(:,n_user)= abs(mean(cap_ma,2));
SE_RO(:,n_user) = abs(mean(cap_ro,2));
SE_FIS(:,n_user)  = abs(mean(cap_fis,2));
end

save('SE_ROMA.mat', 'SE_ROMA');
save('SE_MA.mat', 'SE_MA');
save('SE_RO.mat', 'SE_RO');
save('SE_FIS.mat', 'SE_FIS');
% figure(1);
% p_list = -5:5:30;
% h1=plot(p_list, SE_ROMA ,'-','linewidth',1.5);
% hold on;
% h2=plot(p_list, SE_MA ,'-','linewidth',1.5);
% hold on;
% h3=plot(p_list, SE_RO ,'-','linewidth',1.5);
% hold on;
% h5=plot(p_list, SE_FIS ,'-','linewidth',1.5);
% hold on;
% legend([h1 h2 h3 h5],'ROMA ','MA','RO','FAS ','Interpreter','latex' )
% 
% xlabel('Normalized region size ($\lambda$)','Interpreter','latex')
% ylabel('Spectral Efficiency (bps/Hz)','Interpreter','latex')
% grid on;

