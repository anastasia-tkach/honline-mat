

beta_true = [37.1409 28.436 14.0033 37.0552 20.5967 12.8329 40.4944 23.2687 15.9074 37.9263 23.9032 13.6009 31.9997 19.1319 12.8205 9.72575 3.87206 -7.16046 ...
    25.3963 50.8191 2.28707 8.13206 52.8925 7.5463 -5.93253 49.1046 7.7991 -18.1006 44.8872 3.8978 23.9371 45.629 -17.9329 38.572 7.31059 3.0928 ...
    -10.0083 1.01374 28.5967 41.231 4.05823 48.6988 12 6.73474 1.52305 10 15.6827 10.3118 7.30846 7.0311 8.44143 7.55251 5.85299 5.17427 7.68834 7.64302 ...
    5.67621 5.58432 7.26768 6.97092 5.01217 4.84959 7.84562 6.22559 4.77166 4.18002 9.50548 10.726 10.2172 8.98482 13.2994 13.6244 13.4193 12.9863, 15.4];

num_betas = 75;
num_thetas = 34;
num_iters = 12;


%% Betas
fileID = fopen('E:\Data\sensor-sequences\synthetic\solutions.txt', 'r');
thetas_betas = fscanf(fileID, '%f');
N = length(thetas_betas)/(num_betas + num_thetas);
thetas_betas = reshape(thetas_betas, num_betas + num_thetas, N)';
betas = thetas_betas(:, num_thetas + 1:end);

errors = betas - repmat(beta_true, length(betas), 1);
errors_norm = zeros(size(errors, 1), 1);
for j = 1:length(errors)
    errors_norm(j) = norm(errors(j, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 46:74] + 1));
end
fclose(fileID);

%% Errors
fileID = fopen('E:\Data\sensor-sequences\synthetic\online_weighted_metrics.txt', 'r');
data_errors = fscanf(fileID, '%f');
fclose(fileID);

%% Plot data metric
dark_red = [188, 58, 117]/255;
light_red = 0.5 * [217, 154, 143]/255 + 0.5 * [238, 198, 199]/255;
dark_green = [61, 131, 119]/255;
light_green = [144, 194, 171]/230;

figure_size = [0.25, 0.25, 0.5, 0.6];
figure_borders = [0.07 0.11 0.90 0.84];

for i = length(errors_norm):length(errors_norm)
    f = figure('units', 'normalized', 'outerposition', figure_size); hold on;
    plot((1:i), errors_norm(1:i), 'lineWidth', 2.5, 'color', dark_green);    
    %plot(1:length(data_errors), data_errors, 'lineWidth', 2);
    xlabel('frame number');
    ylabel('\beta - \beta_{true}');
    set(gca,'position', figure_borders, 'units','normalized');
    ylim([0, 42]);
    xlim([0.999, i]);
    set(gca,'fontsize', 15, 'fontname', 'Cambria'); box on;
    
    print(f,['E:\Data\honline-video\Plots\', num2str(i)],'-dpng');
    %close;
end