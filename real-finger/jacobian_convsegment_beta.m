% clc; close all; clear;
% 
% n = 3;
% 
% while (true)
%     
%     while(true)
%         c1 = 0.5 * rand(n ,1);
%         c2 = 0.5 * rand(n ,1);
%         x = (c2 - c1) / norm(c2 - c1);
%         beta = norm(c2 - c1);
%         x1 = rand(1 ,1);
%         x2 = rand(1 ,1);
%         r1 = max(x1, x2);
%         r2 = min(x1, x2);
%         p = rand(n, 1);
%         if norm(c1 - c2) > r1
%             break;
%         end
%         
%     end
%     
%     % we want the point to project on the conic surface, not on the
%     % spherical one
%     [index, q, ~, ~] = projection_convsegment(p, c1, c2, r1, r2, 1, 2);
%     if length(index) == 2, break; end
% end
% 
% p = q;
% display_result({c1, c2}, {p}, {q}, {[1, 2]}, {r1, r2}, true, 0.7, 'small');

function [f, df] = jacobian_convsegment_beta(p, c1, beta, r1, r2, x, variables)

D = length(p);

p_ = p;
c1_ = c1;
beta_ =  beta;
x_ =  x;
r1_ =  r1;
r2_ = r2;

Jnumerical = [];
Janalytical = [];

for var = 1:length(variables)
    variable = variables{var};
    switch variable
        case 'c1', dc1 = eye(D, D); dbeta = zeros(1, D);
            dr1 = zeros(1, D); dr2 = zeros(1, D);
            dp = zeros(D, D); dx = zeros(D, D);
        case 'beta', dc1 = zeros(D, 1); dbeta = ones(1, 1);
            dr1 = zeros(1, 1); dr2 = zeros(1, 1);
            dp = zeros(D, 1); dx = zeros(D, 1);
        case 'r1', dc1 = zeros(D, 1); dbeta = 0;
            dr1 = 1; dr2 = 0; 
            dp = zeros(D, 1); dx = zeros(D, 1);
        case 'r2', dc1 = zeros(D, 1); dbeta = 0;
            dr1 = 0; dr2 = 1;
            dp = zeros(D, 1); dx = zeros(D, 1);
    end
    
    %% c2 = c1 + beta * x
    [u, du] = product_derivative(beta_, dbeta, x_, dx);
    [c2_, dc2] = sum_derivative(c1_, dc1, u, du);
    
    %   O = c2_;
    %   dO = dc2;
    %   O = @(beta) O(c1, beta, r1, r2);
    %   Jnumerical = [Jnumerical, my_gradient(O, beta)]
    %   Janalytical = [Janalytical, dO(c1, beta, r1, r2)]
    
    %% u =  c2 - c1; v =  p - c1;
    %[u, du] = difference_derivative(c2_, dc2, c1_, dc1);
    [v, dv] = difference_derivative(p_, dp, c1_, dc1);
    
    %% t - closest point on the axis, t = c1 + alpha * u;
    [s, ds] = dot_derivative(u, du, v, dv);
    [tn, dtn] = product_derivative(s, ds, u, du);
    [uu, duu] = dot_derivative(u, du, u, du);
    [b, db] = ratio_derivative(tn, dtn, uu, duu);
    [t, dt] =  sum_derivative(c1_, dc1, b, db);
    
    %% omega - lenght of the tangent, omega = sqrt(u' * u - (r1 - r2)^2);
    [r, dr] = difference_derivative(r1_, dr1, r2_, dr2);
    [c, dc] = product_derivative(r, dr, r, dr);
    [omega2, domega2] = difference_derivative(uu, duu, c, dc);
    [omega, domega] = sqrt_derivative(omega2, domega2);
    
    %% delta - size of the correction, % delta =  norm(p - t) * (r1 - r2) / omega;
    [a, da] = difference_derivative(p_, dp, t, dt);
    [b, db] = dot_derivative(a, da, a, da);
    [c, dc] = sqrt_derivative(b, db);
    [deltanum, ddeltanum] = product_derivative(c, dc, r, dr);
    [delta, ddelta] = ratio_derivative(deltanum, ddeltanum, omega, domega);
    
    %% w - correction vector, w = delta * u / norm(u);
    [wnum, dwnum] = product_derivative(delta, ddelta, u, du);
    [unorm, dunorm] = sqrt_derivative(uu, duu);
    [w, dw] = ratio_derivative(wnum, dwnum, unorm, dunorm);
    
    %% s - corrected point on the axis, s = t - w
    [s, ds] =  difference_derivative(t, dt, w, dw);
    
    %% gamma - correction in the direction orthogonal to cone surface, gamma = (r1 - r2) * norm(c2 - t + w)/ norm(u);
    [a, da] = difference_derivative(c2_, dc2, t, dt);
    [b, db] = sum_derivative(a, da, w, dw);
    [c, dc] = dot_derivative(b, db, b, db);
    [gammafactor, dgammafactor] = sqrt_derivative(c, dc);
    [gammanum, dgammanum] =  product_derivative(r, dr, gammafactor, dgammafactor);
    [gamma, dgamma] = ratio_derivative(gammanum, dgammanum, unorm, dunorm);
    
    %% q - the point on the model surface, q = s + (p - s) / norm(p - s) * (gamma + r2);
    
    [a, da] = difference_derivative(p_, dp, s, ds);
    [qfactor, dqfactor] = normalize_derivative(a, da);
    [b, db] = sum_derivative(gamma, dgamma, r2_, dr2);
    [c, dc] = product_derivative(b, db, qfactor, dqfactor);
    [q, dq] = sum_derivative(s, ds, c, dc);
    f = q;

    %% Display result
    switch variable
        case 'c1', df.dc1 = dq;
        case 'beta', df.dbeta = dq;
        case 'r1', df.dr1 = dq;
        case 'r2', df.dr2 = dq;
    end
end










