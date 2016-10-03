function [E1, dE1, ddE1] = simple_problem_egh_1_cumulative(E1_, dE1_, ddE1_, X, Y, T, N, w2)

%% Data term
f1 =  exp(X(N) * T(N))^2 - Y(N);
j1 = zeros(1, N);
j1(N) = 2 * T(N) * exp(X(N) * T(N))^2;
h1 = zeros(N, N);
h1(N, N) = 4 * T(N).^2 * exp(X(N) * T(N)).^2;

E1 = E1_ + f1' * f1;
dE1 = [dE1_, 0] +  2 * f1' * j1;
if N == 1
    ddE1_ = 0;
else
    ddE1_ = [ddE1_, zeros(N - 1, 1); zeros(1, N)];
end
ddE1 = ddE1_ + 2 * j1' * j1 + 2 * f1 * h1;

