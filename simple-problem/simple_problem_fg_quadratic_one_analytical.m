function [F, J, h] = simple_problem_fg_quadratic_one_analytical(x2, x0, y, t, h_, w2, settings)

%% Eliminating x_previous
if ~isempty(h_)
    A = h_(1, 1);
    B = h_(1, 2);
    C = h_(2, 1);
    D = h_(2, 2);
    
    if settings.quadratic_one_marginalization
        W3 = D - C * inv(A) * B;
    else
        W3 = D;
    end
else
    W3 = 0;
end

W2 = w2 * eye(1, 1);

% f1 = @(xx) F1(xx)' * F1(xx);
% o1 = @(xx) f1(xx) + F2(xx)' * F2(xx) + Q(xx)' * Q(xx);
% o2 = @(xx) f1(xx) + w2 * (xx(2) - xx(1))' * (xx(2) - xx(1)) + (xx(1) - x_0)' * W3 * (xx(1) - x_0);
%
% x0 = x_0; x1 = xx(1); x2 = xx(2);
%
% o3 = @(x1, x2) f1(xx) + w2 * (x2 - x1)' * (x2 - x1) + (x1 - x0)' * W3 * (x1 - x0);
%
% o4 = @(x1, x2) f1(xx) + w2 * x2' * x2 - 2 * w2 * x1' * x2 +  w2 * x1' * x1 + x1' * W3 * x1 - 2 * x0' * W3 * x1 + x0' * W3 * x0;
%
% d = @(x1, x2) x0' * W3 * x0 + f1(xx) + w2 * x2' * x2;
%
% o = @(x1, x2) x1' * (W2 + W3) * x1 - 2 * x1' * (W2 * x2 + W3 * x0) + d(x1, x2);
%
% do_dx1 = @(x1, x2) 2 * (W2 + W3) * x1 - 2 * (W2 * x2 + W3 * x0);

x1_opt = @(x2) inv(W2 + W3) * W2 * x2 + inv(W2 + W3) * W3 * x0;
M = inv(W2 + W3) * W2;
m = inv(W2 + W3) * W3 * x0;
x1_opt = @(x2) M * x2 + m;

%disp(do_dx1(x1_opt(x2), x2));

%u = @(x1) o(x1, x2);
%v = my_gradient(u, x1);
%disp([v; do_dx1(x1, x2)]);
%disp([f(xx); o(x1, x2)]);

p = @(x2) f1(xx) + w2 * (x2 - x1_opt(x2))' * (x2 - x1_opt(x2)) + (x1_opt(x2) - x0)' * W3 * (x1_opt(x2) - x0);
I = eye(1, 1);

S1 = @(x2) exp(x2 * t)^2 - y;
dS1_dx2 = @(x2) 2 * t * exp(x2 * t)^2;

S2 = @(x2) sqrt(W2) * (x2 - M * x2 - m);
dS2_dx2 = @(x2) sqrt(W2) * (I - M);

if ~isempty(h_)
    S3 = @(x2) sqrt(W3) * (M * x2 + m - x0);
    dS3_dx2 = @(x2) sqrt(W3) * M;
else
    S3 = @(x2) 0;
    dS3_dx2 = @(x2) 0;
end

S = @(x2) [S3(x2); S1(x2); S2(x2)];
dS = @(x2) [dS3_dx2(x2); dS1_dx2(x2); dS2_dx2(x2)];

s = @(x2) S(x2)' * S(x2);
%disp(s(x2));
%disp([p(x2); s(x2)]);
%disp(' ');

F = S(x2);
J = dS(x2);
h = 2 * dS(x2)' * dS(x2);

h = [1, 1; 1, h];



