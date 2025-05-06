function [Z] = update_z(Z, R, d)

opt        = optimset('Display','off');
% main optimization
z_opt = fmincon(@(Z) objective(Z, R), Z, [], [], [], [], [], [], @(Z) nonlcon(Z,d), opt);
Z     = z_opt;

end


%======function======
% constraint
function [c, ceq] = nonlcon(Z, d)

c = [];
N = size(Z, 1);

for n = 1:N
    for i = [1:n-1, n+1:N]
        c = [c; d^2 - norm(Z(n, :) - Z(i, :))^2];
    end
end

ceq =[];
end

% objective function
function cost = objective(Z, R)

cost = 0;
N    = size(Z, 1);

for n = 1:N
    cost = cost + norm(R(n, :) - Z(n, :))^2;
end

end