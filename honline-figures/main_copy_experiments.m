
clc;
experiment_name = 'last2';

sequence_path = 'E:\Data\sensor-sequences\calibration_new_fingers\';
results_path = 'E:\Data\honline-results\marker_positions_metrics\calibration_new_fingers\';

%% Marker-based metrics
marker_base_metics_forlder = 'markers\';
markers_source_filename = [sequence_path, 'marker_based_metrics'];
markers_target_filename = [results_path, marker_base_metics_forlder];

%copyfile([markers_source_filename, '_batch_offline.txt'], [markers_target_filename, '_batch_offline.txt']);
%copyfile([markers_source_filename, '_batch_online.txt'], [markers_target_filename, '_batch_online.txt']);
%copyfile([markers_source_filename, '_kalman_standard.txt'], [markers_target_filename, '_kalman_standard.txt']);
%copyfile([markers_source_filename, '_kalman_extended.txt'], [markers_target_filename, '_kalman_extended.txt']);
%copyfile([markers_source_filename, '_independent.txt'], [markers_target_filename, '_uncalibrated.txt']);


%% Online metrics
online_metrics_folder = 'online\';
online_source_filename = [sequence_path, 'online_continuous_metrics_'];
online_target_filename = [results_path, online_metrics_folder];

copyfile([online_source_filename, 'batch_offline.txt'], [online_target_filename, 'batch_offline_', experiment_name, '.txt']);
%copyfile([online_source_filename, 'batch_online.txt'], [online_target_filename, 'batch_online_',  experiment_name, '.txt']);
%copyfile([online_source_filename, 'kalman_standard.txt'], [online_target_filename, 'kalman_standard_',  experiment_name, '.txt']);
%copyfile([online_source_filename, 'kalman_extended.txt'], [online_target_filename, 'kalman_extended_',  experiment_name, '.txt']);
%copyfile([online_source_filename, 'independent.txt'], [online_target_filename, 'independent_',  experiment_name, '.txt']);