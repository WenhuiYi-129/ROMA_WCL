function [R_s] = update_rs(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L,Z_S, R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, A)

for m=1:Mh*Mv
    lb(m,1)=-0.5*A;
    lb(m,2)=0;
    lb(m,3)=-0.5*A;
    ub(m,1)=0.5*A;
    ub(m,2)=0;
    ub(m,3)=0.5*A;
end
rho = 5;

% opt        = optimset('Display','off');
opt        = optimoptions('fmincon', 'MaxIter', 100, 'Display','off');
% opt = optimoptions('fmincon', 'MaxIterations', 500);

R_opt = fmincon(@(R_s) objective(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L,Z_S, R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, rho), R_s, [], [], [], [], ...
                               lb, ub, [], opt);
R_s     = R_opt;


end


%======function======
% objective function
function cost = objective(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L,Z_S, R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, rho)
delta=0;
M=Mh*Mv;
N=Nh*Nv;
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
cost = -abs(SE);

for n = 1:M
    cost = cost + rho * (norm(R_s(n, :) - Z_S(n, :))^2);
end

end