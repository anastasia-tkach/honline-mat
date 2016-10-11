%% Display online equivalent
%{
figure_borders = [0.05 0.08 0.93 0.90];
means = zeros(N, B);
vars = zeros(N, B);
trues = zeros(N, B);
%beta_std = zeros(N, B);
for i = 1:length(history)
    means(i, :) = history{i}.betas{i};
    trues(i, :) = beta_true;
    if (run_kalman_filter)
        vars(i, :) = diag(history{i}.P(1:B, 1:B));
    else
        vars(i, :) = diag(history{i}.JtJ((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B));
    end
    %beta_std(i, :) = history{i}.beta_std;
end

figure('units', 'normalized', 'outerposition', [0.1, 0.2, 0.33, 0.7]); hold on;
set(gca,'position', figure_borders, 'units','normalized');
for i = 1:B - 1
    h = subplot(B - 1, 1, i); hold on;
    p = get(h, 'pos'); set(h, 'pos', [0.05, p(2) - 0.06, 0.9, 0.43]);
    set(gca,'XTick',[]);
    plot(1:length(history), means(:, i), '.-', 'lineWidth', 2, 'markersize', 13);
    plot(1:length(history), trues(:, i), 'lineWidth', 2);
    plot(1:length(history), beta_init(i) * ones(length(history), 1), 'lineWidth', 2, 'lineStyle', '-.');
    
    if (run_kalman_filter)
        plot(1:length(history), means(:, i) + vars(:, i).^0.5, 'lineWidth', 2, 'color', [0.55, 0.75, 0.8]);
        plot(1:length(history), means(:, i) - vars(:, i).^0.5, 'lineWidth', 2, 'color', [0.55, 0.75, 0.8]);
    end
    
    %plot(1:length(history), trues(:, i) + 0.5 * beta_variance_threshold * ones(length(history), 1), 'lineWidth', 1, 'color', [0.7, 0.9, 0.5]);
    %plot(1:length(history), trues(:, i) + 0.5 * beta_std(:, i).^2, 'lineWidth', 1, 'color', [0.4, 0.7, 0.5]);
    
    ylim([1.5, 4.5]); xlim([0, length(history)]);
    set(gca, 'fontSize', 13); title(['beta ', num2str(i)]);
end

if exist('video_writer', 'var'), video_writer.close(); end
%}