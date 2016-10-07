function [F, J, H] = test_function(x, D)

m = size(x, 1);
A = reshape(D(1:m * m), m, m);
B = D(m * m + 1:m * m + m)';
c = D(m * m + m + 1);

F = x' * A' * A * x +  B * x + c;
J = 2 * x' * A' * A + B;
H = 2 * A' * A;