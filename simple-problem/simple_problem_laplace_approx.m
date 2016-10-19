function [xx_opt, J, h] = simple_problem_laplace_approx(X0, X_prev, h_, Y, T, N, w2, num_iters)

y_ = Y(N - 1);
t_ = T(N - 1);
y = Y(N);
t = T(N);
x_0 = X_prev(N - 1);
x__0 = X_prev(N - 2);

if N == 3
    [F_, J_, h_] = simple_problem_fg_laplace_approx(X_prev(N - 2: N - 1), [], y_, t_, h_, w2);
end

%% Quadratic approximation
xx0 = X0(N-1:N);

[xx_opt] = my_lsqnonlin(@(xx) simple_problem_fg_laplace_approx(xx, x_0, y, t, h_, w2), xx0, num_iters);
[F, J, h] = simple_problem_fg_laplace_approx(xx_opt, x_0, y, t, h_, w2);

%% Display
% display_posterior_and_laplace_approx(xx_opt, F, J, h, x_0, y, t, h_, w2);

