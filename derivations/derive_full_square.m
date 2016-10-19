clc; clear;
n = 3;
A1 = rand(n, n);
A = A1' * A1;
b = rand(n, 1);
C1 = rand(n, n);
C = C1' * C1;
d = rand(n, 1);
x = rand(n, 1);

f = (A * x + b)' * (A * x + b) + (C * x + d)' * (C * x + d);

f1 = x' * (A' * A + C' * C) * x + 2 * x' * (A * b + C * d) + b' * b + d' * d;

M = sqrtm(A' * A + C' * C);
m = inv(sqrtm(A' * A + C' * C)) * (A * b + C * d);

f2 = x' * M' * M * x + 2 * x' * M * m + m' * m - m' * m + b' * b + d' * d;

k = b' * b + d' * d - m' * m;
f3 = (M * x + m)' * (M * x + m) +  k;

disp([f; f3]);