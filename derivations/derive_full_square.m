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

L = sqrtm(A' * A + C' * C);
l = inv(sqrtm(A' * A + C' * C)) * (A * b + C * d);

f2 = x' * L' * L * x + 2 * x' * L * l + l' * l - l' * l + b' * b + d' * d;

k = b' * b + d' * d - l' * l;
f3 = (L * x + l)' * (L * x + l) +  k;

disp([f; f3]); % f is identical to f3