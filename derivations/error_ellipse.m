clear; clc;
%close all;
%rng default;

% Create some random data
V = randn(2, 2);
sigma = V' * V;
disp(sigma);
mu = randn(2, 1);
N = 5000;
data = mvnrnd(mu, sigma, N);

% Calculate the eigenvectors and eigenvalues
mu = mean(data);
sigma = cov(data);

draw_covariance_matrix(mu', sigma);
%plot(data(:,1), data(:,2), '.', 'color', [0.8, 0.8, 0.8]);