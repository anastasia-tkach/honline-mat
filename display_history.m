function [] = display_history(history, beta_init, beta_variance_threshold, N, B, T)

if isempty(history), return; end

means = zeros(length(history), B + T);
vars = zeros(length(history), B + T);
trues = zeros(length(history), B + T);
beta_std = zeros(length(history), B);
energies = zeros(length(history), B - 1);
for i = 1:length(history)
    means(i, :) = history{i}.mean;
    trues(i, :) = history{i}.true;
    vars(i, :) = diag(history{i}.P);
    beta_std(i, :) = history{i}.beta_std;
    energies(i, :) = history{i}.energy;
end

for i = 1:B - 1
    h = subplot(B - 1, 1, i); hold on;
    p = get(h, 'pos'); set(h, 'pos', [0.05, p(2) - 0.06, 0.9, 0.43]);
    set(gca,'XTick',[]);
    plot(1:length(history), means(:, i), '.-', 'lineWidth', 2, 'markersize', 13);
    plot(1:length(history), trues(:, i), 'lineWidth', 2);
    plot(1:length(history), beta_init(i) * ones(length(history), 1), 'lineWidth', 2, 'lineStyle', '-.');
    
    plot(1:length(history), means(:, i) + vars(:, i).^0.5, 'lineWidth', 2, 'color', [0.55, 0.75, 0.8]);
    plot(1:length(history), means(:, i) - vars(:, i).^0.5, 'lineWidth', 2, 'color', [0.55, 0.75, 0.8]);
    
    %plot(1:length(history), trues(:, i) + 0.5 * beta_variance_threshold * ones(length(history), 1), 'lineWidth', 1, 'color', [0.7, 0.9, 0.5]);
    %plot(1:length(history), trues(:, i) + 0.5 * beta_std(:, i).^2, 'lineWidth', 1, 'color', [0.4, 0.7, 0.5]);
    
    %plot(1:length(history), trues(:, i) - 3 * energies(:, i), 'lineWidth', 1, 'color', [0.8, 0.8, 0.8]);
    
    %legend({'mean', 'true', 'initial', 'mean + std', 'mean - std'}); 
    ylim([1.5, 4.5]); xlim([0, N]);
    set(gca, 'fontSize', 13); title(['beta ', num2str(i)]);
end