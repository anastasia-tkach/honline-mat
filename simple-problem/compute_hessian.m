function [ddE_an] = compute_hessian(X, Y, T, N, w2)

if N == 0
    ddE_an = zeros(0, 0);
    return
end

F1 = @(X, N) exp(diag(T(1:N)) * X(1:N)).^2 - Y(1:N);
dF1 = @(X, N) diag(2 * diag(T(1:N)) * exp(diag(T(1:N)) * X(1:N)).^2);
for i = 1:N, ddF1{i} = @(X, N) diag([zeros(i - 1, 1); 4 * T(i).^2 * exp(diag(T(i)) * X(i)).^2; zeros(N - i, 1)]); end

F2 = @(X, N) X(2:N) - X(1:N-1);
dF2 = @(X, N) [zeros(N - 1, 1), eye(N - 1, N - 1)] - eye(N-1, N);
for i = 1:N, ddF2{i} = @(X, N) zeros(N - 1, N); end

F = @(X, N) [F1(X, N); sqrt(w2) * F2(X, N)];
dF = @(X, N) [dF1(X, N); sqrt(w2) * dF2(X, N)];
E = @(X, N) F(X, N)' * F(X, N);
dE = @(X, N) 2 *  F(X, N)' * dF(X, N);

ddE_an = 2 *  dF(X, N)' * dF(X, N);
for i = 1:N
    ddE_an(i, :) = ddE_an(i, :) + 2 * F(X, N)' * [ddF1{i}(X, N); ddF2{i}(X, N)];
end

