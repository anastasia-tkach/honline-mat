function [E, dE] = simple_problem_fg_for_fminunc(X, Y, T, N, batch_size, w2)

%f = exp(x * t)^2;

%% Data term
F1 = zeros(N, 1);
J1 = zeros(N, N);
for i = max(1, N - batch_size + 1):N
    F1(i) =  exp(X(i) * T(i))^2 - Y(i);
    J1(i, i) = 2 * T(i) * exp(X(i) * T(i))^2;
end

%% Closeness term
F2 = zeros(N - 1, 1);
J2 = zeros(N - 1, N);
for i = max(1, N - batch_size):N - 1
    F2(i) = X(i) - X(i + 1);
    if i > N - batch_size
        J2(i, i) = 1;
    end
    J2(i, i + 1) = -1;
end

F = [F1; sqrt(w2) * F2];
J = [J1; sqrt(w2) * J2];

E = F' * F;
dE = 2 * F' * J;
