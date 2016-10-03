function [ E2, dE2, ddE2] = simple_problem_egh_2_cumulative(E2_, dE2_, ddE2_, X, Y, T, N, w2)

%% Closeness term
f2 = 0;
j2 = zeros(1, N);
h2 = zeros(N, N);
if (N > 1)
    f2 = sqrt(w2) * (X(N - 1) - X(N));
    j2(N - 1) = sqrt(w2) * 1;
    j2(N) = sqrt(w2) * (- 1);
end

E2 = E2_ + f2' * f2;
dE2 = [dE2_, 0] +  2 * f2' * j2;
if N == 1
    ddE2_ = 0;
else
    ddE2_ = [ddE2_, zeros(N - 1, 1); zeros(1, N)];
end
ddE2 = ddE2_ + 2 * j2' * j2 + 2 * f2 * h2;

