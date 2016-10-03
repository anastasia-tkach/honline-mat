function [E1, dE1, ddE1] = simple_problem_fgh_cumulative(E1_, dE1_, ddE1_, X, Y, T, N, w2)

f1 =  exp(X(N) * T(N))^2 - Y(N);
j1 = zeros(1, N);
j1(N) = 2 * T(N) * exp(X(N) * T(N))^2;
h1 = zeros(N, N);
h1(N, N) = 4 * T(N).^2 * exp(X(N) * T(N)).^2;

f2 = 0;
j2 = zeros(1, N);
h2 = zeros(N, N);

if (N > 1)
    f2 = X(N - 1) - X(N);
    j2(N - 1) = 1;
    j2(N) = - 1;
end
f = [f1; sqrt(w2) * f2];
j = [j1; sqrt(w2) * j2];

E1 = E1_ + f' * f;
dE1 = [dE1_, 0] +  2 * f' * j;
if N == 1
    ddE1_ = 0;
else
    ddE1_ = [ddE1_, zeros(N - 1, 1); zeros(1, N)];
end
ddE1 = ddE1_ + 2 * j' * j + 2 * f1 * h1 + 2 * f2 * h2;

