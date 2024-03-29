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


% c1 = centers{1};
% c2 = centers{2};
% r1 = radii{1};
% r2 = radii{2};
% x = rotated_frame{2};

p = q;
display_result({c1, c2}, {p}, {q}, {[1, 2]}, {r1, r2}, false, 0.7, 'small');

D = 3;

arguments = 'c1, beta, r1, r2';
variables = {'beta'};
p_ = @(c1, beta, r1, r2) p;
c1_ = @(c1, beta, r1, r2) c1;
beta_ = @(c1, beta, r1, r2) beta;
x_ = @(c1, beta, r1, r2) x;
r1_ = @(c1, beta, r1, r2) r1;
r2_ = @(c1, beta, r1, r2) r2;

Jnumerical = [];
Janalytical = [];

for var = 1:length(variables)
    variable = variables{var};
    switch variable
        case 'c1', dc1= @(c1, beta, r1, r2) eye(D, D); dbeta = @(c1, beta, r1, r2) zeros(1, D);
            dr1 = @(c1, beta, r1, r2) zeros(1, D); dr2 = @(c1, beta, r1, r2) zeros(1, D);
            dp = @(c1, beta, r1, r2) zeros(D, D);  dx = @(c1, beta, r1, r2) zeros(D, D);
        case 'beta',  dc1= @(c1, beta, r1, r2) zeros(D, 1); dbeta = @(c1, beta, r1, r2) ones(1, 1);
            dr1 = @(c1, beta, r1, r2) zeros(1, 1); dr2 = @(c1, beta, r1, r2) zeros(1, 1);
            dp = @(c1, beta, r1, r2) zeros(D, 1); dx = @(c1, beta, r1, r2) zeros(D, 1);
        case 'r1', dc1= @(c1, beta, r1, r2) zeros(D, 1); dbeta = @(c1, beta, r1, r2) 0;
            dr1 = @(c1, beta, r1, r2) 1; dr2 = @(c1, beta, r1, r2) 0;
            dp = @(c1, beta, r1, r2) zeros(D, 1);  dx = @(c1, beta, r1, r2) zeros(D, 1);
        case 'r2', dc1= @(c1, beta, r1, r2) zeros(D, 1); dbeta = @(c1, beta, r1, r2) 0;
            dr1 = @(c1, beta, r1, r2) 0; dr2 = @(c1, beta, r1, r2) 1;
            dp = @(c1, beta, r1, r2) zeros(D, 1);  dx = @(c1, beta, r1, r2) zeros(D, 1);
    end
    
    %% c2 = c1 + beta * x
    [u, du] = product_handle(beta_, dbeta, x_, dx, arguments);
    [c2_, dc2] = sum_handle(c1_, dc1, u, du, arguments);
    
%     O = c2_;
%     dO = dc2;
%     O = @(beta) O(c1, beta, r1, r2);
%     Jnumerical = [Jnumerical, my_gradient(O, beta)]
%     Janalytical = [Janalytical, dO(c1, beta, r1, r2)]
    
    %% u =  c2 - c1; v =  p - c1;
    %[u, du] = difference_handle(c2_, dc2, c1_, dc1, arguments);
    [v, dv] = difference_handle(p_, dp, c1_, dc1, arguments);
    
    %% t - closest point on the axis, t = c1 + alpha * u;
    [s, ds] = dot_handle(u, du, v, dv, arguments);
    [tn, dtn] = product_handle(s, ds, u, du, arguments);
    [uu, duu] = dot_handle(u, du, u, du, arguments);
    [b, db] = ratio_handle(tn, dtn, uu, duu, arguments);
    [t, dt] =  sum_handle(c1_, dc1, b, db, arguments);
    
    %% omega - lenght of the tangent, omega = sqrt(u' * u - (r1 - r2)^2);
    [r, dr] = difference_handle(r1_, dr1, r2_, dr2, arguments);
    [c, dc] = product_handle(r, dr, r, dr, arguments);
    [omega2, domega2] = difference_handle(uu, duu, c, dc, arguments);
    [omega, domega] = sqrt_handle(omega2, domega2, arguments);
    
    %% delta - size of the correction, % delta =  norm(p - t) * (r1 - r2) / omega;
    [a, da] = difference_handle(p_, dp, t, dt, arguments);
    [b, db] = dot_handle(a, da, a, da, arguments);
    [c, dc] = sqrt_handle(b, db, arguments);
    [deltanum, ddeltanum] = product_handle(c, dc, r, dr, arguments);
    [delta, ddelta] = ratio_handle(deltanum, ddeltanum, omega, domega, arguments);
    
    %% w - correction vector, w = delta * u / norm(u);
    [wnum, dwnum] = product_handle(delta, ddelta, u, du, arguments);
    [unorm, dunorm] = sqrt_handle(uu, duu, arguments);
    [w, dw] = ratio_handle(wnum, dwnum, unorm, dunorm, arguments);
    
    %% s - corrected point on the axis, s = t - w
    [s, ds] =  difference_handle(t, dt, w, dw, arguments);
    
    %% gamma - correction in the direction orthogonal to cone surface, gamma = (r1 - r2) * norm(c2 - t + w)/ norm(u);
    [a, da] = difference_handle(c2_, dc2, t, dt, arguments);
    [b, db] = sum_handle(a, da, w, dw, arguments);
    [c, dc] = dot_handle(b, db, b, db, arguments);
    [gammafactor, dgammafactor] = sqrt_handle(c, dc, arguments);
    [gammanum, dgammanum] =  product_handle(r, dr, gammafactor, dgammafactor, arguments);
    [gamma, dgamma] = ratio_handle(gammanum, dgammanum, unorm, dunorm, arguments);
    
    %% q - the point on the model surface, q = s + (p - s) / norm(p - s) * (gamma + r2);
    
    [a, da] = difference_handle(p_, dp, s, ds, arguments);
    [qfactor, dqfactor] = normalize_handle(a, da, arguments);
    [b, db] = sum_handle(gamma, dgamma, r2_, dr2, arguments);
    [c, dc] = product_handle(b, db, qfactor, dqfactor, arguments);
    [q, dq] = sum_handle(s, ds, c, dc, arguments);
    O = q;
    dO = dq;
    %% Display result
    switch variable
        case 'c1'
            O = @(c1) O(c1, beta, r1, r2);
            Jnumerical = [Jnumerical, my_gradient(O, c1)];
            Janalytical = [Janalytical, dO(c1, beta, r1, r2)];
        case 'beta'
            O = @(beta) O(c1, beta, r1, r2);
            Jnumerical = [Jnumerical, my_gradient(O, beta)];
            Janalytical = [Janalytical, dO(c1, beta, r1, r2)];
        case 'r1'
            O = @(r1) O(c1, beta, r1, r2);
            Jnumerical = [Jnumerical, my_gradient(O, r1)];
            Janalytical = [Janalytical, dO(c1, beta, r1, r2)];
        case 'r2'
            O = @(r2) O(c1, beta, r1, r2);
            Jnumerical = [Jnumerical, my_gradient(O, r2)];
            Janalytical = [Janalytical, dO(c1, beta, r1, r2)];
    end
end
%disp('F = '); disp(O(r2))
Jnumerical
Janalytical
myline(O(beta), O(beta) + Janalytical / norm(Janalytical) * beta, 'r');
myline(O(beta), O(beta) + O(beta) - s(c1,beta,r1,r2), 'b');







