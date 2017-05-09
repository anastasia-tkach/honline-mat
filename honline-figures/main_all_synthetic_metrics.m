beta_true = [37.1409 28.436 14.0033 37.0552 20.5967 12.8329 40.4944 23.2687 15.9074 37.9263 23.9032 13.6009 31.9997 19.1319 12.8205 9.72575 3.87206 ...
    -7.16046 25.3963 50.8191 5.28707 8.13206 52.8925 10.5463 -5.93253 49.1046 10.7991 -18.1006 44.8872 6.8978 23.9371 45.629 -17.9329 38.572 ...
    7.31059 3.0928 -10.0083 1.01374 28.5967 41.231 4.05823 48.6988 12 6.73474 1.52305 15.6827 10.3118 7.30846 7.0311 8.44143 7.55251 ...
    5.85299 5.17427 7.68834 7.64302 5.67621 5.58432 7.26768 6.97092 5.01217 4.84959 7.84562 6.22559 4.77166 4.18002 9.50548 10.726 10.2172 8.98482 13.2994 13.6244 13.4193];

beta_batch = [38.1441 31.1173 12.3628 34.664 21.0186 14.9943 39.0325 23.3399 16.6795 36.0385 23.5144 14.9535 29.2559 19.7434 14.3156 15.6489 4.45096 -6.03128 31.4971 ...
    53.9855 1.29259 14.4743 56.7279 5.90665 0.974826 52.7083 6.48529 -11.0239 47.5729 2.38342 31.4529 47.0504 -10.592 40.6051 8.07634 0.198744 -1.98377...
    0.962445 36.766 50.6443 4.27451 43.8032 13.1798 8.34506 2.16528 15.9981 10.3783 7.20171 7.06595 8.51422 7.61788 5.81748 5.11803 7.76214 7.55892 5.64525...
    5.54565 7.41708 6.92325 4.97765 4.8028 7.99174 6.21545 4.73624 4.14501 10.0618 10.6704 10.0381 9.23011 13.5674 13.7941 13.748];

num_betas = 72;
num_thetas = 0;
num_iters = 6;


%% Data Hmodel
errors = zeros(10, num_betas, length(experiment_names));

for i = 1:length(experiment_names)
    display([data_path, experiment_names{i}, '.txt']);
    fileID = fopen([data_path, experiment_names{i} , '.txt'], 'r');
    
    thetas_betas = fscanf(fileID, '%f');
    N = length(thetas_betas)/(num_betas + num_thetas);
    thetas_betas = reshape(thetas_betas, num_betas + num_thetas, N)';
    %N = 873 * num_iters;
    %N = 3000;
    betas = thetas_betas(start_offset:N, num_thetas + 1:end);
    
    errors(1:size(betas, 1), :, i) = betas - repmat(beta_true, size(betas, 1), 1);
    fclose(fileID);
end;

errors = errors(1:num_iters:end, :, :);

%% Find mean and std
mean_errors = mean(errors, 3);
std_errors = std(errors, [], 3);


%% Plot data metric
figure_size = [0.25, 0.25, 0.5, 0.6];
figure_borders = [0.05 0.08 0.93 0.84];

% [15, 18, 21, 24, 27] + 1 - do not include
% [17, 20, 23, 26, 29] + 1 - do not include
% [30, 32, 34, 36] + 1 - do not include
% [45:71] + 1

beta_indices_cells = {[3, 4, 5, 19], [6, 7, 8, 22], [9, 10, 11, 25], [12, 13, 14, 28]};%, [45:71]};

for f = 1:length(beta_indices_cells)
    figure('units', 'normalized', 'outerposition', figure_size); hold on;
    legend_strings = {};
    i = 0;
    for beta_index = beta_indices_cells{f} + 1
        i = min(i + 1, 7);
        plot(1:N/num_iters, mean_errors(:, beta_index), 'lineWidth', 2, 'color', colors{i});
        %plot(1:N, (beta_batch(beta_index) - beta_true(beta_index)) * ones(N, 1), 'lineWidth', 2, 'color', colors{i});
        legend_strings{end + 1} = num2str(beta_index - 1);
        %plot(1:N, mean_errors(:, beta_index) + std_errors(:, beta_index), 'lineWidth', 1, 'color', colors{i});
        %plot(1:N, mean_errors(:, beta_index) - std_errors(:, beta_index), 'lineWidth', 1, 'color', colors{i});
    end
    %if (f == 1), plot(1:N/num_iters, mean_errors(:, 3 + 1) + mean_errors(:,19 + 1), 'lineWidth', 2, 'color','k'); end
    %if (f == 2), plot(1:N/num_iters, mean_errors(:, 6 + 1) + mean_errors(:,22 + 1), 'lineWidth', 2, 'color','k'); end
    %if (f == 3), plot(1:N/num_iters, mean_errors(:, 9 + 1) + mean_errors(:,25 + 1), 'lineWidth', 2, 'color','k'); end
    %if (f == 4), plot(1:N/num_iters, mean_errors(:, 12 + 1) + mean_errors(:,28 + 1), 'lineWidth', 2, 'color','k'); end
    legend(legend_strings);
    plot(1:N/num_iters, zeros(N/num_iters, 1), 'lineWidth', 1, 'color', 'k');
    set(gca,'position', figure_borders, 'units','normalized');
    ylim([-10, 10]);
    drawnow;
end

