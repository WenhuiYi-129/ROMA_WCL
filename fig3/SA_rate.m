function SA_SE = SA_rate(noise,Nh,Nv,Mh,Mv,lambda,theta_r, phi_r, theta_t,phi_t, path_loss,L,location_u,R_s2, R_r2, M, M2,N, N2,K,cap_temp,cap_before,epsil,alpha,beta,delta)
      

        % 生成所有可能的索引组合
        all_indices_t = nchoosek(1:M2, M);
        all_indices_r = nchoosek(1:N2, N);
fre = 200;
         SA_SE=-1;
        % 遍历所有组合
        
for i = 1:size(all_indices_t, 1)
            selected_indices_t = all_indices_t(i, :); % 当前组合的索引
            % 提取选择的坐标点
           R_s = R_s2(selected_indices_t,:);
           R_r=zeros(N,3,K);
        j=1;
while (j<= fre)
    j=j+1;
    b=toc;
fprintf('SA beginning, SA_SE = %f, times is %.4f\n',SA_SE,b);
    for u=1:K
        t=randi([1,size(all_indices_r,1)]);
        selected_indices_r = all_indices_r(t, :); % 当前组合的索引
       R_r(:,:,u)=R_r2(selected_indices_r,:,u);
    end
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
    W_precoding(:,:,u) = pinv(B)*H_channel(:,:,u);
    W_precoding(:,:,u) =W_precoding(:,:,u)/norm(W_precoding(:,:,u), 'fro');
end
SE= SE_compute(H_channel, W_precoding, K, N, noise);
cap_temp1 = SE;
if cap_temp1>SA_SE
    SA_SE =cap_temp1;
end
end
end

end