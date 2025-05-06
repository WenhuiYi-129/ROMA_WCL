function r_matrix = position(x0, y0, z0, alpha, beta,Mh,Mv,d,RS)
%%ZS表示天线单元在面阵上的相对位置，即不考虑旋转角度所在的坐标
M = Mh*Mv;
r_matrix = zeros(M,3);

if RS==0
for m=1:M
    m1 = mod(m-1, Mh);
    m2 = floor((m-1)/Mh);
    delta_x = (m1-(Mh-1)/2)*d*cos(alpha)-(m2-(Mv-1)/2)*d*sin(beta)*sin(alpha);
    delta_y = (m1-(Mh-1)/2)*d*sin(alpha)+(m2-(Mv-1)/2)*d*sin(beta)*cos(alpha);
    delta_z = (m2-(Mv-1)/2)*d*cos(beta);
    r_matrix(m,1) = x0+delta_x;
    r_matrix(m,2) = y0+delta_y;
    r_matrix(m,3) = z0+delta_z;
end
else
for m=1:M
    delta_x = RS(m,1)*cos(alpha)-RS(m,3)*sin(beta)*sin(alpha);
    delta_y =  RS(m,1)*sin(alpha)+RS(m,3)*sin(beta)*cos(alpha);
    delta_z = RS(m,3)*cos(beta);
    r_matrix(m,1) = x0+delta_x;
    r_matrix(m,2) = y0+delta_y;
    r_matrix(m,3) = z0+delta_z;
end
end

end