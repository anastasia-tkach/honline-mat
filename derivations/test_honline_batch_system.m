clc; clear; close all;

path = 'C:/Users/tkach/Desktop/Test/';

num_thetas = 29;
num_betas = 30;
num_thetas_latent = 4;
num_betas_latent = 1;
num_poses = 3;

num_parameteres = num_thetas + num_betas + num_thetas_latent + num_betas_latent;
num_batch_parameters = num_poses * (num_thetas + num_thetas_latent) + num_betas + num_betas_latent;

system_lhs = cell(num_poses, 1);
system_rhs = cell(num_poses, 1);

for i = 1:num_poses
    % read system-lhs
    fileID = fopen([path, 'system-', num2str(i - 1), '-lhs.txt'], 'r');
    system_lhs{i} = fscanf(fileID, '%f');
    system_lhs{i} = reshape(system_lhs{i}, num_parameteres, num_parameteres)';
    
    % read system-rhs
    fileID = fopen([path, 'system-', num2str(i - 1), '-rhs.txt'], 'r');
    system_rhs{i} = fscanf(fileID, '%f');    
    
    figure; spy(system_lhs{i});
end

% read system-lhs
fileID = fopen([path, 'batch-system-lhs.txt'], 'r');
batch_system_lhs = fscanf(fileID, '%f');
batch_system_lhs = reshape(batch_system_lhs, num_batch_parameters, num_batch_parameters)';

% read system-lhs
fileID = fopen([path, 'batch-system-rhs.txt'], 'r');
batch_system_rhs = fscanf(fileID, '%f');

figure; spy(batch_system_lhs);

return;

%% Thetas system

thetas_system_lhs = [system_lhs(1:num_thetas, 1:num_thetas), ...
    system_lhs(1:num_thetas, num_betas + num_thetas + 1:num_betas + num_thetas + num_thetas_latent);
    system_lhs(num_thetas + num_betas + 1:num_thetas + num_betas + num_thetas_latent, 1:num_thetas),...
    system_lhs(num_thetas + num_betas + 1:num_thetas + num_betas + num_thetas_latent, num_thetas + num_betas + 1:num_thetas + num_betas + num_thetas_latent)];

thetas_system_rhs = [system_rhs(1:num_thetas); system_rhs(num_thetas + num_betas + 1:num_thetas + num_betas + num_thetas_latent)];

batch_thetas_system_lhs = batch_system_lhs(1:num_thetas + num_thetas_latent, 1:num_thetas + num_thetas_latent);
batch_thetas_system_rhs = batch_system_rhs(1:num_thetas + num_thetas_latent);

thetas_solution = thetas_system_lhs \ thetas_system_rhs;
batch_thetas_solution = batch_thetas_system_lhs \ batch_thetas_system_rhs;

%disp([thetas_solution, batch_thetas_solution]);

%% Betas system

betas_system_lhs = [system_lhs(num_thetas + 1:num_thetas + num_betas, num_thetas + 1:num_thetas + num_betas), ...
    system_lhs(num_thetas + 1:num_thetas + num_betas, num_thetas + num_betas + num_thetas_latent + 1:end);
    system_lhs(num_thetas + num_betas + num_thetas_latent + 1:end, num_thetas + 1:num_thetas + num_betas),...
    system_lhs(num_thetas + num_betas + num_thetas_latent + 1:end, num_thetas + num_betas + num_thetas_latent + 1:end)];

betas_system_rhs = [system_rhs(num_thetas + 1:num_thetas + num_betas); system_rhs(num_thetas + num_betas + num_thetas_latent + 1:end)];

batch_betas_system_lhs = batch_system_lhs(num_thetas + num_thetas_latent + 1:end, num_thetas + num_thetas_latent + 1:end);
batch_betas_system_rhs = batch_system_rhs(num_thetas + num_thetas_latent + 1:end);

betas_solution = betas_system_lhs \ betas_system_rhs;
batch_betas_solution = batch_betas_system_lhs \ batch_betas_system_rhs;

%disp([betas_solution, batch_betas_solution]);

%return

%% Solve the systems

solution = system_lhs \ system_rhs;

batch_solution = batch_system_lhs \ batch_system_rhs;

solution = solution([1:num_thetas, num_thetas + num_betas + 1:num_thetas + num_betas + num_thetas_latent, ...
    num_thetas + 1:num_thetas + num_betas, end - num_betas_latent + 1:end]);

disp([solution, batch_solution]);

