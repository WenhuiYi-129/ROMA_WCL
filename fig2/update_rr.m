function [R_rr] = update_rr(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L,Z_Rr, R_s, R_rr, u, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, A)

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

R_opt = fmincon(@(R_rr) objective(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L,Z_Rr, R_s, R_rr, u, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, rho), R_rr, [], [], [], [], ...
                               lb, ub, [], opt);
R_rr    = R_opt;


end


%======function======
% objective function
function cost = objective(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L,Z_Rr, R_s, R_rr, u, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda, rho)
delta=0;
M=Mh*Mv;
N=Nh*Nv;
H_channel = zeros(M,N,K);
W_precoding = zeros(M,N,K);
B=zeros(M,M);
for i =1:K
    r_t = position(0, 0, 0, alpha(1), beta(1),Mh,Mv,delta,R_s);
    if i~=u
    r_r = position(location_u(i,1),location_u(i,2),location_u(i,3), alpha(i+1),beta(i+1),Nh,Nv,delta,R_r(:,:,i));
    else
       r_r = position(location_u(i,1),location_u(i,2),location_u(i,3), alpha(i+1),beta(i+1),Nh,Nv,delta,R_rr);
    end  
    H_channel(:,:,i)=Channel(L, path_loss(i,:), theta_t, phi_t, theta_r, phi_r,...
 lambda, r_t, r_r, M, N);
    B=B+H_channel(:,:,i)*H_channel(:,:,i)';
end
for i=1:K
    W_precoding(:,:,i) = inv(B)*H_channel(:,:,i);
    W_precoding(:,:,u) =W_precoding(:,:,u)*(sqrt(p_max)/norm(W_precoding(:,:,u), 'fro'));
end
SE= SE_compute(H_channel, W_precoding, K, N, noise);
cost = -abs(SE);

for n = 1:N
    cost = cost + rho * (norm(R_rr(n, :) - Z_Rr(n, :))^2);
end

end