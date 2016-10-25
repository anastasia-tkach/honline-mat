clc; clear;
options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt');


n = 2;
m = 1;

x = rand(n, 1);

Q = ones(1, n);
y = 1/n * Q * x;
z = 1/n * (Q' * Q) * x;

%% Version 1
F1 = @(x) x - 1/n * (Q' * Q) * x;
x1 = lsqnonlin(F1, x, [], [], options);

%% Version 2
F2 = @(xy) xy(1:n) -  Q' * xy(n + 1:end);
xy = [x; y];
xy2 = lsqnonlin(F2, xy, [], [], options);

%% Version 3
h_sqrt = (eye(n, n) -  1/n * (Q' * Q));
F3 = @(x) h_sqrt * x;
x3 = lsqnonlin(F3, x, [], [], options);

mu = [0; 0];
epsilon = 0.001;
h_sqrt = h_sqrt + sqrt(epsilon) * eye(n, n);
sigma = inv(h_sqrt' * h_sqrt);
draw_covariance_matrix(mu, sigma, 1);
title(['\epsilon = ', num2str(epsilon)]);

% options = optimoptions(@fminunc,'Algorithm','quasi-newton');
% [x_opt, fval, exitflag, output] = fminunc(E2_opt4, x, options);

disp(x1);
disp(xy2(1:n));
disp(x3);

% disp(x);
% disp(mean(x));
% disp(y);
% disp(z);



