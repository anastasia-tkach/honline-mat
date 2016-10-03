rng(46546);

n = 5;
m = 2;

J = rand(n, m);
F = randn(n, 1);
r = 0.25;
R = diag(r * ones(n, 1));

RJ = R^(-0.5) * J;
RF = R^(-0.5) * F;
I = diag([1, 2]);

x1 = (J' * inv(R) * J + I) \ (J' * inv(R) * F);
x2 = (J' * J + r * I) \ (J' * F);
x3 = (RJ' * RJ + I) \ (RJ' * RF);
x4 = (J' * J) \ (J' * F);
disp([x1, x2, x3, x4]);

% when r is bigger, damping is bigger. r = 1 equavalent to damping = 1

