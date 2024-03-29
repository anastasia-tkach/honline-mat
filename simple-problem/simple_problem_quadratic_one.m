function [xx_opt, J, h] = simple_problem_quadratic_one(X0, X_prev, h_, Y, T, N, w2, settings)

y_ = Y(N - 1);
t_ = T(N - 1);
y = Y(N);
t = T(N);
x_0 = X_prev(N - 1);

if N == 3
    [~, ~, h_] = simple_problem_fg_quadratic_one(X_prev(N - 1), X_prev(N - 2), y_, t_, [], w2);   
end

%% Quadratic approximation
%xx0 = X0(N-1:N);
x2 = X0(N);
[x2_opt] = my_lsqnonlin(@(x2) simple_problem_fg_quadratic_one(x2, x_0, y, t, h_, w2), x2, settings.num_iters);

W2 = w2 * eye(1, 1);

A = h_(1, 1); B = h_(1, 2); C = h_(2, 1); D = h_(2, 2);
W3 = D;

x1_opt = inv(W2 + W3) * W2 * x2_opt + inv(W2 + W3) * W3 * x_0;

xx_opt = [x1_opt; x2_opt];

[F, J, h] = simple_problem_fg_quadratic_one(x2_opt, x_0, y, t, h_, w2);




