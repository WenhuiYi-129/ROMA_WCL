clc
clear 


% para
L=15; % number of channel paths
f_c = 2.1e9;
c = 3e8;
lambda = c/f_c;               % \lambda
A_max  = lambda*3;           % max of A
A_min  = lambda*1;           % min of A
N_A    = 5;                 % number of different A

K = 4;  %user_num
Nh      = 3;     
Nv      = 3;  
N = Nh*Nv;
% N
Mh      = 3; 
Mv     =3;
M = Mh*Mv;
% M
%Noise figure
  noisefigure=9;
  %communication bandwidth
  B=20e6;
  noisepowerdbm=-174+10*log10(B)+noisefigure;%%noise power
  noise=db2pow(noisepowerdbm)/1000;
% SNR        = 15;             % SNR dB
% SNR_linear = 10^(SNR/10);    % SNR
% noise = 1/SNR;
p_max_list = [-10 0 10];

epsil      = 10^-3; 
d          = lambda/2; 
Nter_list      =1:10:301;
% cap        = zeros(N_A, Nter);     % capacity

%%用户分布范围
x_max = 5;
x_min = -5;
y_max = 5;
y_min = -5;
z_max = 5;
z_min = -5;
tic; 
location_u = zeros(K,3);
%%用户位置
location_u(:,1)=(x_max-x_min)*rand(K,1)+x_min;
location_u(:,2)=(y_max-y_min)*rand(K,1)+y_min;
location_u(:,3)=(z_max-z_min)*rand(K,1)+z_min;

%%计算路径损耗
path = zeros(K,1);
path_loss = zeros(K,L);

PATH=0;
for u=1:K
distances = sqrt(location_u(u,1)^2+location_u(u,2)^2+location_u(u,3));
[path(u),~] =  functionChannelgain(distances);
path(u) = sqrt(db2pow(path(u)));
PATH = PATH+path(u)^2*M*N;
z = (randn(1, L) + 1i * randn(1, L));  % 随机生成复数数列，实部和虚部均为正态分布
% 计算当前数列的2范数
current_norm = norm(z, 'fro');

path_loss(u,:) = path(u)*z/current_norm; % 生成复数序列，模为1
end


%%需要优化的参数有天线相对平面的位置，两个旋转角度
ENDpoint=zeros(length(Nter_list),length(p_max_list));
ENDpoint_DE = zeros(length(Nter_list),length(p_max_list));
ENDMAX = zeros(length(Nter_list),1);
ENDH = zeros(length(Nter_list),1);
theta_t = rand(L,1)*pi;                                          % \theta_t
phi_t   = rand(L,1)*pi;                                          % \phi_t
theta_r = rand(L,1)*pi;                                          % \theta_r
phi_r   = rand(L,1)*pi;                                          % \phi_r
n_A=4;
A = (A_max-A_min)/(N_A-1)*(n_A-1) + A_min;
delta =A/2;%%天线间距初始 
% Initialization of alpha and beta for M=N=4, d=lambda/2
alpha = rand(K+1,1)*pi;
beta = rand(K+1,1)*pi;
for i_p = 1:length(p_max_list)
p_max = p_max_list(i_p);
p_max = db2pow(p_max)/1000;



R_s = position(0, 0, 0, 0,0,Mh,Mv,delta,0); % shape 2*4
R_r = zeros(N,3,K);
for u=1:K
R_r(:,:,u) = position(location_u(u,1), location_u(u,2), location_u(u,3), 0,0,Nh,Nv,delta,0);
end
% Initialization of Z
Z_S = (2*rand(M, 3)-1)*A;
Z_R = (2*rand(N,3,K)-1)*A;
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
    W_precoding(:,:,u) =W_precoding(:,:,u)*sqrt(p_max)/norm(W_precoding(:,:,u), 'fro');
end
SE= SE_compute(H_channel, W_precoding, K, N, noise);
% SE_base(n_A,iter)=SE;
cap_temp = SE;
cap_before = -1;

NP=400;
x = zeros(M+K*N+2*K+2, 3,NP);
flag=0;
for list = 1:length(Nter_list)
    Nter = Nter_list(list);
    SE_list = zeros(Nter,1);
    [SE_max,x]=DE4D(NP, Nter, x, p_max,noise,...
    location_u,Mh,Mv,Nh,Nv,K,L,theta_r, phi_r, ...
    theta_t,phi_t,path_loss, lambda, A,flag);
    flag=1;
for iter = 1:Nter
    b=toc;
cap_before = cap_temp;


% update {r_s}, fix {r_r}, alpha, beta and {z_n}
[R_s] = update_rs(p_max, noise,location_u,Mh,Mv,Nh,Nv,K,L,Z_S, R_s, R_r, theta_r, phi_r, theta_t,phi_t, alpha,beta,path_loss, lambda, A);
% update {z_n}, fix Q and {r_n}
[Z_S] = update_z(Z_S, R_s, d);
% update {r_r}, fix {r_s}, alpha, beta and {z_n}
% for u=1:K
% R_rr = R_r(:,:,u);
% Z_Rr = Z_R(:,:,u);
% R_rr = update_rr(p_max, noise,location_u,Mh,Mv,Nh,Nv,K,L,Z_Rr, R_s, R_rr, u, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, A);
% R_r(:,:,u) = R_rr;
% [Z_Rr ] = update_z(Z_Rr , R_rr , d);
% Z_R(:,:,u) = Z_Rr;
% end

[alpha] = update_alpha(noise,location_u,Mh,Mv,Nh,Nv,K,L,R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, A);


[beta] = update_beta(noise,location_u,Mh,Mv,Nh,Nv,K,L, R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, A);


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
    W_precoding(:,:,u) =W_precoding(:,:,u)*sqrt(p_max)/norm(W_precoding(:,:,u), 'fro');
end
SE= SE_compute(H_channel, W_precoding, K, N, noise);
cap_temp = SE;
fprintf('now iter is %d/%d/301, pieces %d, times is %.4f\n', iter,Nter, cap_temp, b);
% HGAIN=0;
% for u=1:K
% HGAIN = HGAIN+norm(H_channel(:,:,u),'fro')^2;
% end
end
ENDpoint_DE(list,i_p)=abs(SE_max);
ENDpoint(list,i_p)=abs(cap_temp);
end
end
% SE_average = abs(mean(cap,2));
% SE_base = abs(mean(SE_base,2));
figure();
h1=plot(Nter_list, ENDpoint(:,1));
hold on;
h2=plot(Nter_list, ENDpoint(:,2));
hold on;
h3=plot(Nter_list, ENDpoint(:,3));
hold on ;
h5=plot(Nter_list, ENDpoint_DE(:,1));
hold on ;
h6=plot(Nter_list, ENDpoint_DE(:,2));
hold on ;
h7=plot(Nter_list, ENDpoint_DE(:,3));
hold on ;
% h1=plot(((A_max-A_min)/(N_A-1)*(0:N_A-1) + A_min)/lambda, SE_average,'-','linewidth',1.5);
% hold on;
% h2=plot(((A_max-A_min)/(N_A-1)*(0:N_A-1) + A_min)/lambda, SE_base,'-','linewidth',1.5);
% hold on;
legend('$p$=-10 dBm','$p$=0 dBm','$p$=10 dBm','Interpreter','latex' )

xlabel('Number of iterations','Interpreter','latex')
ylabel('Spectral Efficiency (bps/Hz)','Interpreter','latex')
grid on;
