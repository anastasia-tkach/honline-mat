function [F, J] = simple_problem_fg_quadratic_all(x, X_, Y, T, E1, dE1, ddE1, N, w2)

K = [zeros(1, N - 2), 1];

%X_opt = M1 * x + M2;

M1 = - (2 * w2 * K' * K  +  ddE1 + ddE1') \ (- 2 * w2 * K)';
M2 = - (2 * w2 * K' * K  +  ddE1 + ddE1') \ (dE1 - 2 * X_' * ddE1)';

m3 = dE1 * M1 + 2 * M1' * ddE1 * M2 - 2 * X_' * ddE1 * M1;
m4 = - 2 * X_' * ddE1 *  M2 + E1 - dE1 * X_ + X_' * ddE1 * X_ + dE1 * M2 + M2' * ddE1 * M2;

[U, S, ~] = eig(ddE1); S1 = S.^0.5;
if any(imag(S1)), disp('imaginary'); end
M5 = S1' * U' * M1;
if (M5 == 0), m6 = 0;
else m6 = (-0.5 * m3 / M5)'; end

E2_lin = @(x) norm(exp(T(N) * x).^2 - Y(N))^2 + ...
    w2 * norm(x - [zeros(1, N - 2), 1] * (M1 * x + M2))^2 + ...
    x' * (M1' * U * S1 * S1' * U' * M1) * x + m3 * x  + m4;

F_ = @(x) [exp(T(N) * x).^2 - Y(N); sqrt(w2) * (x - [zeros(1, N - 2), 1] * (M1 * x + M2)); M5 * x - m6];

J_ = @(x) [2 * T(N) * exp(T(N) * x).^2; sqrt(w2) * (1 - [zeros(1, N - 2), 1] * M1); M5];

%v = my_gradient(F_, x); disp([J_(x)'; v']);

F = F_(x);
J = J_(x);



