function [F, J, h] = simple_problem_fg_quadratic_one(x2, x0, y, t, h_, w2, settings)

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

M = inv(W2 + W3) * W2;
m = inv(W2 + W3) * W3 * x0;
x1_opt = M * x2 + m;

I = eye(1, 1);

F1 = exp(x2 * t)^2 - y;
dF1_dx2 = 2 * t * exp(x2 * t)^2;

F2 = sqrtm(W2) * (x2 - M * x2 - m);
dF2_dx2 = sqrtm(W2) * (I - M);

if ~isempty(h_)
    F3 = sqrtm(W3) * (M * x2 + m - x0);
    dF3_dx2 = sqrtm(W3) * M;
else
    F3 = 0;
    dF3_dx2 = 0;
end

F = [F3; F1; F2];
J = [dF3_dx2; dF1_dx2; dF2_dx2];

h = 2 * J' * J;

h = [1, 1; 1, h];



