function [capavg]= AS10(Nter,A_max,A_min,N_A,lambda)
% clear;
% clc
% Nter=10;                % Monte Carlo simulations迭代次数
Lt = 10;                    % L_t发射天线角度数量
Lr = 10;                    % L_r
SNR = 15;                   % SNR dB
% lambda = 0.05;              % \lambda
% A_max = lambda*3;           % max of A
% A_min = lambda*1;           % min of A
% N_A = 11;                   % number of different A
N = 4;                      % N
M = 4;                      % M
SNR_linear = 10^(SNR/10);   % SNR
PAPC = (SNR_linear/N)*ones(N,1);
epsil = 10^-3;              % \epsilon_1
eps = 10^-6;             % \epsilon_2
d = lambda/2;               % D

cap = zeros(N_A, Nter);     % capacity
%% 

for iter = 1:Nter
    theta_t = rand(Lt,1)*pi;                                        % \theta_t
    phi_t = rand(Lt,1)*pi;                                          % \phi_t
    angle_t = [sin(theta_t).*cos(phi_t) cos(theta_t)];
    theta_r = rand(Lr,1)*pi;                                        % \theta_r
    phi_r = rand(Lr,1)*pi;                                          % \phi_r
    angle_r = [sin(theta_r).*cos(phi_r) cos(theta_r)]; 
    Sigma = diag(1/sqrt(2*Lr) * (randn(Lr,1) + 1j*randn(Lr,1)));    % \Sigma
%% 

    for n_A =1:N_A
        A = A_max;                  % lambda-3lambda 改变可移动天线的移动区域大小

        % Initialization of r and t for M=N=4, d=lambda/2
        rx = [-A/4 -A/4 A/4 A/4];                                   % x_r
        ry = [-A/4 A/4 A/4 -A/4];                                   % y_r
        tx = [-A/4 -A/4 A/4 A/4];                                   % x_t
        ty = [-A/4 A/4 A/4 -A/4];                                   % y_t
        
        H = exp(1j*2*pi*angle_r*[rx;ry]/lambda)'*Sigma*exp(1j*2*pi*angle_t*[tx;ty]/lambda);  % H
        [cap_temp,Q_half] = singleantenna(H, PAPC, eps);
        W = H*Q_half;                                               %  
       %% 
       
        
        cap_before = -1;
        %opt = optimset('Display','off','TolCon',1e-12);
                               
        rx = [-(3*A)/4 -(3*A)/4 -A/4 -A/4 A/4 A/4 (3*A)/4 (3*A)/4]; % x_r
        ry = [-A/4 A/4 -A/4 A/4 A/4 -A/4 A/4 -A/4]; % y_r
        tx = [-(3*A)/4 -(3*A)/4 -A/4 -A/4 A/4 A/4 (3*A)/4 (3*A)/4];                                  % x_t
        ty = [-A/4 A/4 -A/4 A/4 A/4 -A/4 A/4 -A/4];

        % 生成所有可能的索引组合
        all_indices = nchoosek(1:length(rx), 4);

        % 遍历所有组合
        for i = 1:size(all_indices, 1)
            selected_indices = all_indices(i, :); % 当前组合的索引

            % 提取选择的坐标点
            selected_rx = rx(selected_indices);
            selected_ry = ry(selected_indices);
            for j = 1:size(all_indices, 1)
                selected_indices = all_indices(j, :); % 当前组合的索引

                % 提取选择的坐标点
                selected_tx = tx(selected_indices);
                selected_ty = ty(selected_indices);
                H = exp(1j*2*pi*angle_r*[selected_rx; selected_ry]/lambda)'*Sigma*exp(1j*2*pi*angle_t*[selected_tx;selected_ty]/lambda);
                [cap_temp,Q_half] = singleantenna(H, PAPC, eps);
                W = H*Q_half;
                if cap_temp>cap_before
                    cap_before = cap_temp;
                end
            end
        end
        cap(n_A,iter) = cap_before;
    end
end

capavg = mean(cap,2);%返回每一种A的均值end
% figure;
% plot(capavg);
%capavg1=capavg(1:10:110);
end
