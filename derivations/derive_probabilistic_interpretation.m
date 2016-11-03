rng default;
clc; clear; %close all;
options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'display','off');
R = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];

chisquare_val = 2.4477;
n = 2;

sigma1 = diag([0.1, 2]);
%sigma1 = diag([2, 2]);
mu1 = zeros(n, 1);

sigma2 = R(-pi/4) * diag([2, 0.1]) * R(-pi/4)';
mu2 = 0.5 * randn(2, 1);

num_runs = 200;

%% optimize
X0 = mvnrnd(mu1, sigma1, num_runs);
X1 = zeros(num_runs, n);
X2 = zeros(num_runs, n);

w0 = 0;
w1 = 1;
w2 = 1;
w3 = 1;

for run_index = 1:num_runs
    
    x = X0(run_index, :)';   
    x0 = x;
    xx = [x; x];
    
    F = @(xx)  [sqrt(w0) * (xx(1:2) - x0); ...
        sqrt(w1) * inv(sqrtm(sigma1)) * (xx(1:2) -  mu1); 
        sqrt(w2) * (xx(1:2) - xx(3:4)); ...
        sqrt(w3) * inv(sqrtm(sigma2)) * (xx(3:4) -  mu2)];
    xx = lsqnonlin(F, xx, [], [], options);
    
    X1(run_index, :) = xx(1:2)';
    X2(run_index, :) = xx(3:4)';
    
end

%% display

figure; hold on; axis equal;

%% first distirbution
[ellipse_points] = get_covarince_elipse(sigma1, chisquare_val);
plot(ellipse_points(:,1) + mu1(1), ellipse_points(:,2) + mu1(2), '-', 'lineWidth', 2, 'color', [1.0 0.45 0.3]);
scatter(X0(:, 1), X0(:, 2), 20, [206, 173, 209]/255, 'o', 'filled');

%% second distribution
[ellipse_points] = get_covarince_elipse(sigma2, chisquare_val);
plot(ellipse_points(:,1) + mu2(1), ellipse_points(:,2) + mu2(2), '-', 'lineWidth', 2, 'color', [34, 177, 76]/255);
scatter(X1(:, 1), X1(:, 2), 20, [1.0 0.45 0.3], 'o', 'filled');
scatter(X2(:, 1), X2(:, 2), 20, [34, 177, 76]/255, 'o', 'filled');

%% connecting distribution
sigma12 = inv(sqrtm(w2 * eye(n, n)));
[ellipse_points] = get_covarince_elipse(sigma12, chisquare_val);
plot(ellipse_points(:,1) + X1(1, 1), ellipse_points(:,2) + X1(1, 2), '-', 'lineWidth', 2, 'color', [0.8, 0.8, 0.8]);

 
 