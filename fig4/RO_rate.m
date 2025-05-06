function RO_SE = RO_rate(p_max,d,epsil,cap_temp, cap_before,iter,Nter, n_A, noise,...
    location_u,Mh,Mv,Nh,Nv,K,L, R_s, R_r, theta_r, phi_r, ...
    theta_t,phi_t, alpha,beta,path_loss, lambda, A)
M=Mh*Mv;
N=Nh*Nv;

while(cap_temp-cap_before > epsil)
% for i=1:100
cap_before = cap_temp;

[alpha] = update_alpha(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L,R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, A);


[beta] = update_beta(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L,R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, A);


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
cap_temp = SE;
b=toc;
fprintf('now iter is %d/%d, space is %d, pieces %d/%d, times is %.4f\n', iter,Nter, n_A, cap_temp, cap_before, b);

end
RO_SE = cap_before ;

end