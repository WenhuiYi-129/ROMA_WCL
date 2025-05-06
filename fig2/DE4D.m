function [SE_max,x]=DE4D(NP, G, x, p_max,noise,...
    location_u,Mh,Mv,Nh,Nv,K,L,theta_r, phi_r, ...
    theta_t,phi_t,path_loss, lambda, A,flag)
% NP=35;%%种群数目
D=1+K+2;%%变量的维数
% G=300;%%最大进化代数
F0=0.5;%%初始变异算子
CR=0.2;%%交叉算子

M=Mh*Mv;
N=Nh*Nv;
%%需要优化的变量有，基站M个天线相对位置（M×3),K个用户N个天线相对位置（N×3),
%%1个基站加上U个用户的旋转角度，alpha(K+1×1），beta(K+1×1）
%%所以定义种群x的每个个体为（(M+KN+2K+2)×3）的矩阵
%%因此定义种群x为四维矩阵

%%天线相对位置的上限和下限
Xs(1,:)=[0.5*A,0,0.5*A];
Xx(1,:)=[-0.5*A,0,-0.5*A];
%%旋转角度alpha的上限和下限
Xs(2,:)=[pi, 0, 0];
Xx(2,:)=[0, 0, 0];
%%旋转角度beta的上限和下限
Xs(3,:)=[pi,0,0];
Xx(3,:)=[0,0,0];
if flag ==0 %%如果flag不等于0则继承上次迭代的结果作为初值
%%赋初值%%
x=zeros(M+K*N+2*K+2, 3,NP);%%初始种群
v=zeros(M+K*N+2*K+2, 3,NP);%%变异种群
u=zeros(M+K*N+2*K+2, 3,NP);%%选择种群
x(1:M,:,:)=rand(M,3,NP).*(Xs(1,:)-Xx(1,:))+Xx(1,:); 
x((M+1):(M+K*N),:,:)=rand(K*N,3,NP).*(Xs(1,:)-Xx(1,:))+Xx(1,:); 
x((M+K*N+1):(M+K*N+K+1),1,:)=rand(K+1,1,NP).*(Xs(2,1)-Xx(2,1))+Xx(2,1); 
x((M+K*N+K+2):(M+K*N+2*K+2),1,:)=rand(K+1,1,NP).*(Xs(3,1)-Xx(3,1))+Xx(3,1); 
end
%%计算目标函数%%
for m=1:NP
    R_s = x(1:M,:,m);
    R_r = zeros(N,3,K);
    for k=1:K
    R_r(:,:,k)=x((M+(k-1)*N+1):(M+k*N),:,m);
    end
    alpha=x((M+K*N+1):(M+K*N+K+1),1,m);
    beta=x((M+K*N+K+2):(M+K*N+2*K+2),1,m);
   SE = SE_general(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L, R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda);
   SE1(m)=SE;
end
[SE_max,index]=max(SE1);
%%跟踪原始种群中的最优个体

for gen=1:G
b=toc;
fprintf('now iter is %d/%d, pieces %d, times is %.4f\n', gen,G, SE_max, b);

  %%%%%变异操作%%%%%%
  %%%%%自适应变异算子%%%%%
  lamb=exp(1-G/(G+1-gen));%自适应差分进化算法
  F=F0*2^(lamb);%动态的设置变异算子
  %%%%%r1,r2,r3和m 互不相同
  for m=1:NP
      r1=randi([1,NP],1,1);
      while (r1==m)
          r1=randi([1,NP],1,1);
      end
      r2=randi([1,NP],1,1);
      while (r2==m)||(r2==r1)
          r2=randi([1,NP],1,1);
      end
      r3=randi([1,NP],1,1);
      while(r3==m)||(r3==r1)||(r3==r2)
          r3=randi([1,NP],1,1);
      end
      v(:,:,m)=x(:,:,r1)+F*(x(:,:,r2)-x(:,:,r3));%%计算变异个体
  end
  %%%%%%交叉操作%%%%%%
  r=randi([1,D],1,1);  %取r防止u==x
  cr=rand(1);
  if (cr<=CR)||(r==1)
      u(1:M,:,:)=v(1:M,:,:);
  else
      u(1:M,:,:)=x(1:M,:,:);
  end
 for k=1:K
      cr = rand(1);
 if (cr<=CR)||(r==1+k)
      u((M+(k-1)*N+1):(M+k*N),:,:,:)=v((M+(k-1)*N+1):(M+k*N),:,:,:);
 else
      u((M+(k-1)*N+1):(M+k*N),:,:,:)=x((M+(k-1)*N+1):(M+k*N),:,:,:);
 end
 end
  cr=rand(1);
  if (cr<=CR)||(r==2+K)
      u((M+K*N+1):(M+K*N+K+1),1,:)=v((M+K*N+1):(M+K*N+K+1),1,:);
  else
      u((M+K*N+1):(M+K*N+K+1),1,:)=x((M+K*N+1):(M+K*N+K+1),1,:);
  end
    cr=rand(1);
  if (cr<=CR)||(r==3+K)
      u((M+K*N+K+2):(M+K*N+2*K+2),1,:)=v((M+K*N+K+2):(M+K*N+2*K+2),1,:);
  else
      u((M+K*N+K+2):(M+K*N+2*K+2),1,:)=x((M+K*N+K+2):(M+K*N+2*K+2),1,:);
  end
  %%%%%%边界条件的处理%%%%%%
  for n=1:M+K*N+2*K+2
      for m=1:NP
         if (n<=M+K*N)
              if (u(n,:,m)<Xx(1,:))|(u(n,:,m)>Xs(1,:)) 
                   u(n,:,m)=rand*(Xs(1,:)-Xx(1,:))+Xx(1,:);
              end
         else 
               if (u(n,:,m)<Xx(2,:))|(u(n,:,m)>Xs(2,:)) 
                   u(n,:,m)=rand*(Xs(2,:)-Xx(2,:))+Xx(2,:);
              end
         end

      end
  end
  %%%%%%选择操作%%%%%
for m=1:NP
    R_s = u(1:M,:,m);
    R_r = zeros(N,3,K);
    for k=1:K
    R_r(:,:,k)=u((M+(k-1)*N+1):(M+k*N),:,m);
    end
    alpha=u((M+K*N+1):(M+K*N+K+1),1,m);
    beta=u((M+K*N+K+2):(M+K*N+2*K+2),1,m);
   SE = SE_general(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L, R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda);
   SE2(m)=SE;
end
for m=1:NP
    if   SE2(m)>  SE1(m)
        x(:,:,m)=u(:,:,m);
    end
end
for m=1:NP
    R_s = x(1:M,:,m);
    R_r = zeros(N,3,K);
    for k=1:K
    R_r(:,:,k)=x((M+(k-1)*N+1):(M+k*N),:,m);
    end
    alpha=x((M+K*N+1):(M+K*N+K+1),1,m);
    beta=x((M+K*N+K+2):(M+K*N+2*K+2),1,m);
   SE = SE_general(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L, R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda);
   SE1(m)=SE;
end
[SE_max,index]=max(SE1);
% %%跟踪原始种群中的最优个体
% for D_i=1:D
% trace(gen+1,D_i)=x(D_i,index);
% end
end
% [~,index]=max(rank1);
alpha1_end=x(1,index);
beta1_end=x(2,index);
alpha2_end=x(3,index);
beta2_end=x(4,index);
disp("4Dend");
end


%% computer SE
function SE = SE_general(p_max,noise,location_u,Mh,Mv,Nh,Nv,K,L, R_s, R_r, theta_r, phi_r, theta_t,phi_t,alpha,beta,path_loss, lambda)
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
    W_precoding(:,:,u) =W_precoding(:,:,u)*sqrt(p_max)/norm(W_precoding(:,:,u), 'fro');
end
SE= SE_compute(H_channel, W_precoding, K, N, noise);

end
