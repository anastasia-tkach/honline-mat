function [xx_opt, J, h] = simple_problem_laplace_approx(X0, X_prev, h_, Y, T, N, w2, num_iters)

y_ = Y(N - 1);
t_ = T(N - 1);
y = Y(N);
t = T(N);
x_0 = X_prev(N - 1);
x__0 = X_prev(N - 2);

%{
f1 =  exp(x_0 * t_)^2 - y_;
j1 = 2 * t_ * exp(x_0 * t_)^2;
h1 = 4 * t_.^2 * exp(x_0 * t_).^2;

f2 = sqrt(w2) * (x_0 - x__0);
j2 = sqrt(w2) * 1;
h2 = 0;

E1_ = f1' * f1; E2_ = f2' * f2;
dE1_ = 2 * f1' * j1; dE2_ = 2 * f2' * j2;
ddE1_ = 2 * j1' * j1 + 2 * f1 * h1; 
ddE2_ = 2 * j2' * j2 + 2 * f2 * h2;
E_ = E1_ + E2_; dE_ = dE1_ + dE2_; ddE_ = ddE1_ + ddE2_;
h_ = [0, 0; 0, ddE_];
%}

%%{
if N == 3
    [F_, J_, h_] = simple_problem_fg_laplace_approx(X_prev(N - 2: N - 1), [], y_, t_, h_, w2);
end
%%}

%% Quadratic approximation
xx0 = X0(N-1:N);

[xx_opt] = my_lsqnonlin(@(xx) simple_problem_fg_laplace_approx(xx, x_0, y, t, h_, w2), xx0, num_iters);
[F, J, h] = simple_problem_fg_laplace_approx(xx_opt, x_0, y, t, h_, w2);


