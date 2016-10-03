function [F, J] = simple_problem_fg_last_n(X, X_prev, history, y, t, N, w2, w3, m)


%% Data term
F1 = exp(X(N) * t)^2 - y;
J1 = zeros(1, N);
J1(N) = 2 * t * exp(X(N) * t)^2;

%% Closeness term
F2 = zeros(N - 1, 1);
J2 = zeros(N - 1, N);
for i = max(1, N - m + 1):N - 1
    F2(i) = sqrt(w2) * (X(i) - X(i + 1));
    J2(i, i) = sqrt(w2);
    J2(i, i + 1) = - sqrt(w2);
end

%% History term
F3 = zeros(N - 1, 1);
J3 = zeros(N - 1, N);
for i = max(1, N - m + 1):N - 1
    A = sqrt(history{i}.JtJ(i, i));
    F3(i) = sqrt(w3) * A * (X(i) - X_prev(i));
    J3(i, i) = sqrt(w3) * A;
end

%% All terms
F = [F1; F2; F3];
J = [J1; J2; J3];

%disp([xx(1), xx(2), x_prev, F'])




