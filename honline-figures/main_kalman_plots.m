
%% Middle middle
%{
sequence_name = 'calib_7_u0.65_w1.1';
beta_index = 7;
y_shift = -2;
y_max_limit = 20;
%}

%% Ring bottom
%{
sequence_name = 'calib_9_u0.7_w1.3';
beta_index = 9;
y_shift = 1;
y_max_limit = 16;
%}

%% Middle z
%{
sequence_name = 'calib_23_u0.65_w1.1';
beta_index = 23;
y_shift = 0;
y_max_limit = 12;
%}

%% Thumb middle
%{
sequence_name = 'calib1_u1.08_w1.1';
beta_index = 1;
y_shift = -2;
y_max_limit = 27;
%}

%% Ring y
%%{
sequence_name = 'calib25_u0.7_w1.3';
beta_index = 25;
y_shift = 0;
y_max_limit = 12;
%%}


%% Beta true

betas_true = [37.1409 28.436 14.0316 37.0552 20.5967 12.8331 40.4944 23.2687 15.9073 37.9263 23.9032 13.6011 31.9997,... 
    19.1319 12.8207 9.72575 3.87206 -7.16046 25.3963 50.8191 2.28707 8.13206 52.8925 7.5463 -5.93253, ...
49.1046 7.7991 -18.1006 44.8872 3.8978 23.9371 45.629 -17.9329 38.572 7.31059 3.0928 -10.0083 1.01374, ...
28.5967 41.231 4.02408 48.6212 12 6.7374 1.52305 10.0241 15.6827 10.3118 7.30846 7.0311 8.44143 7.55251, ...
5.85299 5.17427 7.68834 7.64302 5.67621 5.58432 7.26768 6.97092 5.01217 4.84959 7.84562 6.22559, ...
4.77166 4.18002 9.50548 10.726 10.2172 8.98482 13.2994 13.6244 13.4193 12.9863 15.4];

betas_true(2) = 27.736;

num_betas = 75;
beta_index = beta_index + 1;
beta_true = betas_true(beta_index);

data_root = 'E:\Data\sensor-sequences\calib\';
data_path = [data_root, sequence_name, '\'];
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

certainties_scaling_factor = y_max_limit/2;

f = figure('units', 'normalized', 'outerposition', figure_size); hold on;
for i = 2:1:length(measured_values)
    clf; hold on;
    plot((1:i), zeros(i, 1), 'lineWidth', 3, 'color', [0.75, 0.75, 0.75], 'linestyle', '-');
    
    yyaxis left
    plot_values = measured_values(1:i) - beta_true * ones(i, 1);
    plot((1:i), plot_values, 'lineWidth', 3, 'color', dark_red, 'lineStyle', '-');   
    plot_values = estimated_values(1:i) - beta_true * ones(i, 1);
    plot((1:i), plot_values, 'lineWidth', 3, 'color', dark_green, 'lineStyle', '-');    
    ylim([-3 * y_max_limit/4 - y_shift, y_max_limit/4 - y_shift]);
    set(gca, 'ycolor', 0.85 * dark_red); ylabel('\beta - \beta_{gt} mm');
    
    yyaxis right   
    plot_certainties = estimated_certainties(1:i) ./ max(estimated_certainties(1:i) / max(measured_certainties(1:i)));
    plot((1:i), plot_certainties, 'lineWidth', 3, 'color', light_green, 'lineStyle', '-');  
    plot_certainties = measured_certainties(1:i);    
    plot((1:i), plot_certainties, 'lineWidth', 3, 'color', light_red, 'lineStyle', '-');   
    ylim([0, 2 * max(measured_certainties(1:i))]); ylabel('\partial \partial F \\ \partial \partial \beta');
    set(gca, 'ycolor', 0.85 * light_red); 
    
    xlabel('frame number');    
    set(gca,'position', figure_borders, 'units','normalized');   
    xlim([1, i]);
    set(gca,'fontsize', 15, 'fontname', 'Cambria'); 
    box on; set(gca,'linewidth', 2); 
    set(gca, 'xcolor', [0.3, 0.3, 0.3]);
    %drawnow;
    print(f,['E:\Data\honline-video\FRAMES\calib\', sequence_name, '\plot\', num2str(i)],'-dpng');
    %close;
end