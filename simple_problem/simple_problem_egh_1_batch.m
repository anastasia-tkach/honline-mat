function [E1, dE1, ddE1] = simple_problem_egh_1_batch(X, Y, T, N, batch_size, w2)

%% Data term
F1 = zeros(N, 1);
dF1 = zeros(N, N);
for i = max(1, N - batch_size + 1):N
    F1(i) =  exp(X(i) * T(i))^2 - Y(i);
    dF1(i, i) = 2 * T(i) * exp(X(i) * T(i))^2;
end
for i = 1:N, ddF1{i} = diag([zeros(i - 1, 1); 4 * T(i).^2 * exp(diag(T(i)) * X(i)).^2; zeros(N - i, 1)]); end

E1 = F1' * F1;
dE1 = 2 * F1' * dF1;

if N == 0, ddE1 = zeros(0, 0); return; end

ddE1 = 2 *  dF1' * dF1;
for i = 1:N
    ddE1(i, :) = ddE1(i, :) + 2 * F1' * ddF1{i};
end






