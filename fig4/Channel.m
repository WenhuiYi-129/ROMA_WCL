function H_channel=Channel(path_num, path_loss, theta_s, phi_s, theta_r, phi_r,...
    lambda, r_t, r_r, M, N)

H_channel = zeros(M,N);
A_s = zeros(M,1);
A_r = zeros(1,N);
for l=1:path_num
k_s = (2*pi/lambda)*[cos(theta_s(l))*cos(phi_s(l)), cos(theta_s(l))*sin(phi_s(l)), sin(theta_s(l))];
k_r = (2*pi/lambda)*[cos(theta_r(l))*cos(phi_r(l)), cos(theta_r(l))*sin(phi_r(l)), sin(theta_r(l))];
for m=1:M
    r_tm = r_t(m,:);
    A_s(m,1)=exp(1j*dot(k_s,r_tm));
end

for n=1:N
    r_rn = r_r(n,:);
    A_r(1,n) = exp(-1j*dot(k_r,r_rn));
end


    H_channel=H_channel+sqrt(1/path_num)*path_loss(l)*A_s*A_r;
end

end