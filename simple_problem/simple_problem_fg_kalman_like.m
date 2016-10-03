function [F, J] = simple_problem_fg_kalman_like(x, x_prev, y, t, JtJ, N, w2)

%% Data term
F1 =  exp(x * t)^2 - y;
J1 = 2 * t * exp(x * t)^2;

%% Closeness term
F2 = 0;
J2 = 0;
if (N > 1)
    F2 = sqrt(w2) * sqrt(JtJ) * (x - x_prev);
    J2 = sqrt(w2) * sqrt(JtJ) * 1;
end

%% All terms
F = [F1; F2];
J = [J1; J2];


