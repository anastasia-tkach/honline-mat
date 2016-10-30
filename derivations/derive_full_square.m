clc; clear;
n = 3;
A1 = rand(n, n);
A = A1' * A1;
b = rand(n, 1);
C1 = rand(n, n);
C = C1' * C1;
d = rand(n, 1);
x = rand(n, 1);

%% Derivation
w2 = 2;
x0 = rand(n, 1);
W2 = w2 * eye(n, n);

J = rand(n, n);
W3 = J' * J;

M = inv(W2 + W3) * W2;
m = inv(W2 + W3) * W3 * x0;

I = eye(n, n);

A = sqrtm(W2) * (I - M);
b = sqrt(W2) * m;
C = sqrtm(W3) * M;
d = sqrtm(W3) * (m - x0);

%% Full square

f = (A * x + b)' * (A * x + b) + (C * x + d)' * (C * x + d);

f1 = x' * (A' * A + C' * C) * x + 2 * x' * (A * b + C * d) + b' * b + d' * d;

L = sqrtm(A' * A + C' * C);
l = inv(sqrtm(A' * A + C' * C)) * (A * b + C * d);


f2 = x' * L' * L * x + 2 * x' * L * l + l' * l - l' * l + b' * b + d' * d;

k = b' * b + d' * d - l' * l;
f3 = (L * x + l)' * (L * x + l) +  k;

%disp([f; f3]); % f is identical to f3



o = inv(A' * A + C' * C) * (A * b + C * d);

f4 = (L * (x + o))' * (L * (x + o)) +  k;
f5 = (x + o)' * (L' * L) * (x + o) +  k;
f6 = (x + o)' * (A' * A + C' * C) * (x + o) +  k;

f7 = (x + o)' *  ((I - M)' * W2 * (I - M) + M' * W3 * M) * (x + o) +  k;
f8 = (x + o)' *  ( W2 - 2 * w2^2 * inv(W2 + W3)' + w2^2 * inv(W2 + W3)' * (W2 + W3) * inv(W2 + W3) ) * (x + o) +  k;

f9 = (x + o)' *  ( w2 * I - w2^2 * inv(W2 + W3)) * (x + o) +  k;


disp([f; f9]); 


