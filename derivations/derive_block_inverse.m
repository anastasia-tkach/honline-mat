close all; clc; clear;
m = 3;
J = randn(2 * m, 2 * m);
A = J' * J;
B = inv(A);

a = A(1:m, 1:m);
b = A(1:m, m + 1:2 * m);
c = A(m + 1:2 * m, 1:m);
d = A(m + 1:2 * m, m + 1:2 * m);

%figure; imagesc(A * B - eye(2 * m, 2 * m)); colorbar;

s11 = inv(a - b * inv(d) * c);
s12 = -inv(a - b * inv(d) * c) * b * inv(d);
s21 = -inv(d - c * inv(a) * b) * c * inv(a);
s22 = inv(d - c * inv(a) * b);

S = [s11, s12; s21, s22];

%figure; imagesc(A * S - eye(2 * m, 2 * m)); colorbar;

r = s22 - s21 * inv(s11) * s12;

disp(r);
disp(inv(d));

