function SE= SE_compute(H_total, W_total, user_num, N, noise)
SE = 0;
for u =1:user_num
    H_channel = H_total(:,:,u)';
    E=zeros(N,N);
    for i = 1:user_num
        W_i = W_total(:,:,i);
        B_i = H_channel*W_i;
        if i==u
            B_u = H_channel*W_i;
        end
        E = E+B_i*B_i';
    end
    E = E-B_u*B_u'+noise*eye(N);
    SE_u=log2(det(eye(N)+B_u'*inv(E)*B_u));
    SE = SE+SE_u;
end
SE=SE/user_num;
end