%sequence_name = 'stability';
sequence_name = 'parameterwise';

input_path = 'saved-variables/probabilistic-interpretation/';

load([input_path, sequence_name, '/results_history']);
load([input_path, sequence_name, '/covariance_history']);
load([input_path, sequence_name, '/results_history_batch']);
load([input_path, sequence_name, '/results_history_batch_independent']);
load([input_path, sequence_name, '/results_history_uniform']);

options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'display','off');
R = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
chisquare_val = 2.4477;
Q = ones(1, 2); W4 = (eye(2, 2) -  1/2 * (Q' * Q));
n = 2;
w0 = 1;
w1 = 1;
w2 = 1;
w3 = 1;
w4 = 1;
zero_frame_index = 1;
first_frame_index = 2;
second_frame_index = 3;

figure; hold on; axis equal;

%% set up parameters
X0 = squeeze(results_history(:, zero_frame_index, 1:2));
X1 = zeros(num_runs, n);
X2 = zeros(num_runs, n);
%X_test = squeeze(results_history_batch(:, second_frame_index, 1:2));
%X_test = squeeze(results_history_batch_independent(:, second_frame_index, 1:2));
X_test = squeeze(results_history_uniform(:, second_frame_index, 1:2));

invalid_indices = [];
for run_index = 1:settings.num_runs
    disp(run_index);
    sigma1 = squeeze(covariance_history(run_index, first_frame_index, 1:2, 1:2));
    mu1 = squeeze(results_history(run_index, first_frame_index, 1:2));
    
    sigma2 = squeeze(covariance_history(run_index, second_frame_index, 1:2, 1:2));
    mu2 = squeeze(results_history(run_index, second_frame_index, 1:2));
    
    if any(isinf(sigma1(:))) || any(isinf(sigma2(:)))
        invalid_indices = [invalid_indices; run_index];
        continue;
    end
    
    x = mu1;
    x0 = X0(run_index, :)';
    xx = [x; x];
    
    F = @(xx)  [sqrt(w0) * (xx(1:2) - x0); ...
        sqrt(w1) * inv(sqrtm(sigma1)) * (xx(1:2) -  mu1); ...
        sqrt(w2) * (xx(1:2) - xx(3:4)); ...
        sqrt(w3) * inv(sqrtm(sigma2)) * (xx(3:4) -  mu2); ...
        sqrt(w4) * W4 * xx(3:4);
    ];
    xx = lsqnonlin(F, xx, [], [], options);
    
    X1(run_index, :) = xx(1:2)';
    X2(run_index, :) = xx(3:4)';
    
    %% display
    
    % fist distribution
    [ellipse_points] = get_covarince_elipse(sigma1, chisquare_val);
    plot(ellipse_points(:,1) + mu1(1), ellipse_points(:,2) + mu1(2), '-', 'lineWidth', 2, 'color', [255, 215, 196]/255);
    
    % second distribution
    [ellipse_points] = get_covarince_elipse(sigma2, chisquare_val);
    plot(ellipse_points(:,1) + mu2(1), ellipse_points(:,2) + mu2(2), '-', 'lineWidth', 2, 'color', [185, 215, 174]/255);
    
    % connecting distribution
    sigma12 = inv(sqrtm(w2 * eye(n, n)));
    [ellipse_points] = get_covarince_elipse(sigma12, chisquare_val);
    plot(ellipse_points(:,1) + X1(1, 1), ellipse_points(:,2) + X1(1, 2), '-', 'lineWidth', 2, 'color', [0.8, 0.8, 0.8]);
end

%% remove invalind inidices

X0(invalid_indices, :) = [];
X1(invalid_indices, :) = [];
X2(invalid_indices, :) = [];
X_test(invalid_indices, :) = [];

%% data points
%scatter(X0(:, 1), X0(:, 2), 20, [206, 173, 209]/255, 'o', 'filled');
%scatter(X1(:, 1), X1(:, 2), 20, [1.0 0.45 0.3], 'o', 'filled');
scatter(X2(:, 1), X2(:, 2), 20, [34, 177, 76]/255, 'o', 'filled');

%% plot X_test
scatter(X_test(:, 1), X_test(:, 2), 20, 'r', 'o', 'filled');


xlim([-1, 7]); ylim([-1, 7]);