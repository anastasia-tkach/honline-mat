function [E2, dE2, ddE2] = simple_problem_egh_2_batch(X, Y, T, N, batch_size, w2)

%% Closeness term
F2 = zeros(N - 1, 1);
dF2 = zeros(N - 1, N);
for i = max(1, N - batch_size):N - 1
    F2(i) = sqrt(w2) * (X(i) - X(i + 1));
    if i > N - batch_size
        dF2(i, i) = sqrt(w2) * 1;
    end
    dF2(i, i + 1) = sqrt(w2) * (-1);
end
for i = 1:N, ddF2{i} = zeros(N - 1, N); end

E2 = F2' * F2;
dE2 = 2 * F2' * dF2;

if N == 0, ddE2 = zeros(0, 0); return; end

ddE2 = 2 *  dF2' * dF2;
for i = 1:N
    ddE2(i, :) = ddE2(i, :) + 2 * F2' * ddF2{i};
end




