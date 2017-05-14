
sequence_name = 'filippe\';

data_path = ['E:\Data\sensor-sequences\online_dataset\', sequence_name];
num_betas = 75;

fileID = fopen([data_path, 'estimated_certainties_new.txt'], 'r');
estimated_certainties = fscanf(fileID, '%f'); fclose(fileID);
N = length(estimated_certainties)/(num_betas);
estimated_certainties = reshape(estimated_certainties, num_betas, N)';

fileID = fopen([data_path, 'measured_certainties_new.txt'], 'r');
measured_certainties = fscanf(fileID, '%f'); fclose(fileID);
N = length(measured_certainties)/(num_betas);
measured_certainties = reshape(measured_certainties, num_betas, N)';

% fileID = fopen([data_path, 'estimated_values.txt'], 'r');
% estimated_values = fscanf(fileID, '%f'); fclose(fileID);
% N = length(estimated_values)/(num_betas);
% estimated_values = reshape(estimated_values, num_betas, N)';
% 
% fileID = fopen([data_path, 'measured_values.txt'], 'r');
% measured_values = fscanf(fileID, '%f'); fclose(fileID);
% N = length(measured_values)/(num_betas);
% measured_values = reshape(measured_values, num_betas, N)';

%% Spesify relevant parameters

betas_thumb = [16, 0, 1, 2] + 1;
max_values_thumb = [7, 6, 3, 0.5];

betas_index = [19, 3, 4, 5] + 1;
max_values_index = [9, 3.5, 2.5, 0.3];

betas_middle = [22, 6, 7, 8] + 1;
max_values_middle = [9, 4, 2.5, 0.5];

betas_ring = [25, 9, 10, 11] + 1;
max_values_ring = [7, 4, 2.5, 0.5];

betas_pinky = [28, 12, 13, 14] + 1;
max_values_pinky = [6, 2, 1.5, 0.4];

beta_indices = [betas_thumb, betas_index, betas_middle, betas_ring, betas_pinky];
max_values = [max_values_thumb, max_values_index, max_values_middle, max_values_ring, max_values_pinky];


%% Display
dark_red = [205, 105, 136]/255;
light_red = 0.5 * [217, 154, 143]/255 + 0.5 * [238, 198, 199]/255;
dark_green = [114, 186, 169]/255;
light_green = [144, 194, 171]/220;

figure_size = [0.25, 0.25, 0.5, 0.6];
figure_borders = [0.07 0.11 0.85 0.84];

y_max_limit = 25;
y_shift = 2;
certainties_scaling_factor = y_max_limit/2;

%[2, 121, 323, 460, 621];

%% Measured
%%{
y_measured_limit = 2.5;
f = figure('units', 'normalized', 'outerposition', figure_size); hold on;
for i = 559:1:length(measured_certainties)
    clf;
    plot((1:i), zeros(i, 1), 'lineWidth', 3, 'color', [0.75, 0.75, 0.75]);
    
    measured_certainties_frame = measured_certainties(i, beta_indices);   
    measured_certainties_frame = measured_certainties_frame./max_values;
    bar(measured_certainties_frame, 'facecolor', dark_red, 'edgecolor', 'none');
    
    %xlabel('\beta index');
    set(gca,'position', figure_borders, 'units','normalized'); axis off;
    %set(gca,'fontsize', 15, 'fontname', 'Cambria'); set(gca, 'ycolor', [1, 1, 1]);
    %set(gca,'linewidth', 4); set(gca, 'xcolor', [0.7, 0.7, 0.7]); box on; 
    xlim([0, length(beta_indices) + 1]); ylim([0, y_measured_limit]); set(gca,'ytick',[]); set(gca,'xtick',[]); set(gca, 'ycolor', [1, 1, 1])
    
    drawnow; %pause(0.05);
    
    print(f,['E:\Data\honline-video\FRAMES\online_dataset\', sequence_name, 'h-measured\', num2str(i)],'-dpng');
end
close;
%%}
%% Estimated
%{
f = figure('units', 'normalized', 'outerposition', figure_size); hold on;
for i = 1:1:length(estimated_certainties)
    clf;
    plot((1:i), zeros(i, 1), 'lineWidth', 3.5, 'color', [0.75, 0.75, 0.75]);
    
    estimated_certainties_frame = estimated_certainties(i, beta_indices);   
    estimated_certainties_frame = estimated_certainties_frame./max_values;
   
    bar(estimated_certainties_frame, 'facecolor', dark_green, 'edgecolor', 'none');    
   
    set(gca,'position', figure_borders, 'units','normalized');
    set(gca,'fontsize', 15, 'fontname', 'Cambria'); box on;
    
    xlim([0, length(beta_indices) + 1]); ylim([0, 2.2]); axis off;   set(gca,'color', 'w');
    set(gca,'ytick',[]); set(gca,'xtick',[]); set(gca, 'ycolor', [1, 1, 1]);
    set(gca,'linewidth', 4); drawnow;
    
    print(f,['E:\Data\honline-video\FRAMES\online_dataset\', sequence_name, 'h-estimated\', num2str(i)],'-dpng');
end
disp(min(estimated_certainties_frame));
close
%}