function [F, J] = simple_problem_fg_last_two(xx, x_prev, y, t, JtJ, N, w2, w3)

x = xx(1);
x_ = xx(2);

%% Data term
F1 = exp(x * t)^2 - y;
J1 = [2 * t * exp(x * t)^2, 0];

%% Closeness term
F2 = 0;
J2 = [0, 0];
if (N > 1)
    F2 = sqrt(w2) * (x - x_);
    J2 = sqrt(w2) * [1, - 1];
end

%% History term

F3 = 0;
J3 = [0, 0];
if (N > 1)
    F3 = sqrt(w3) * sqrt(JtJ) * (x_ - x_prev);
    J3 = sqrt(w3) * [0, sqrt(JtJ)];
end

%% All terms
F = [F1; F2; F3];
J = [J1; J2; J3];

%disp([xx(1), xx(2), x_prev, F'])




