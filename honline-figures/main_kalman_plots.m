
data_path = 'E:\Data\sensor-sequences\kalman\calib_9_u0.7_w1.3\';
beta_true = 37.92;%40.5;
beta_index = 9 + 1;
num_betas = 75;

fileID = fopen([data_path, 'estimated_certainties.txt'], 'r');
estimated_certainties = fscanf(fileID, '%f'); fclose(fileID);
N = length(estimated_certainties)/(num_betas);
estimated_certainties = reshape(estimated_certainties, num_betas, N)';
estimated_certainties = estimated_certainties(:, beta_index);

fileID = fopen([data_path, 'measured_certainties.txt'], 'r');
measured_certainties = fscanf(fileID, '%f'); fclose(fileID);
N = length(measured_certainties)/(num_betas);
measured_certainties = reshape(measured_certainties, num_betas, N)';
measured_certainties = measured_certainties(:, beta_index);

fileID = fopen([data_path, 'estimated_values.txt'], 'r');
estimated_values = fscanf(fileID, '%f'); fclose(fileID);
N = length(estimated_values)/(num_betas);
estimated_values = reshape(estimated_values, num_betas, N)';
estimated_values = estimated_values(:, beta_index);

fileID = fopen([data_path, 'measured_values.txt'], 'r'); 
measured_values = fscanf(fileID, '%f'); fclose(fileID);
N = length(measured_values)/(num_betas);
measured_values = reshape(measured_values, num_betas, N)';
measured_values = measured_values(:, beta_index);

%% Display
dark_red = [188, 58, 117]/255;
light_red = 0.5 * [217, 154, 143]/255 + 0.5 * [238, 198, 199]/255;
dark_green = [61, 131, 119]/255;
light_green = [144, 194, 171]/230;

figure_size = [0.25, 0.25, 0.5, 0.6];
figure_borders = [0.07 0.11 0.85 0.84];

y_max_limit = 25;
y_shift = 2;
certainties_scaling_factor = y_max_limit/2;

for i = length(measured_values):length(measured_values)
    f = figure('units', 'normalized', 'outerposition', figure_size); hold on;
    
    plot((1:i), zeros(i, 1), 'lineWidth', 3, 'color', [0.75, 0.75, 0.75]);
    
    yyaxis left
    plot_values = measured_values(1:i) - beta_true * ones(i, 1);
    plot((1:i), plot_values, 'lineWidth', 3, 'color', dark_red, 'lineStyle', '-');   
    plot_values = estimated_values(1:i) - beta_true * ones(i, 1);
    plot((1:i), plot_values, 'lineWidth', 3, 'color', dark_green, 'lineStyle', '-');    
    ylim([-3 * y_max_limit/4 - y_shift, y_max_limit/4 - y_shift]);
    set(gca, 'ycolor', [0.2, 0.2, 0.2]); ylabel('\beta');
    
    yyaxis right   
    plot_certainties = estimated_certainties(1:i) ./ max(estimated_certainties(1:i) / max(measured_certainties(1:i)));
    plot((1:i), plot_certainties, 'lineWidth', 3, 'color', light_green, 'lineStyle', '-');  
    plot_certainties = measured_certainties(1:i);    
    plot((1:i), plot_certainties, 'lineWidth', 3, 'color', light_red, 'lineStyle', '-');   
    ylim([0, 2 * max(measured_certainties(1:i))]); ylabel('\partial \partial F \\ \partial \partial \beta');
    set(gca, 'ycolor', [0.2, 0.2, 0.2]); 
    
    xlabel('frame number');    
    set(gca,'position', figure_borders, 'units','normalized');   
    xlim([2, i]);
    set(gca,'fontsize', 15, 'fontname', 'Cambria'); box on;
    box on; set(gca,'linewidth', 1.5); 
    set(gca, 'xcolor', [0.2, 0.2, 0.2]);
    
    print(f,['E:\Data\honline-video\Plots\', num2str(i)],'-dpng');
    %close;
end