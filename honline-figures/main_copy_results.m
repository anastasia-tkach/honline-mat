

experiment_name = 'batch_offline_reload';

sequence_path = 'E:\Data\sensor-sequences\calibration_new_fingers\';
results_path = 'E:\Data\honline-results\marker_positions_metrics\calibration_new_fingers\';

%% Marker-based metrics
%marker_base_metics_forlder = 'markers\';
%markers_source_filename = [sequence_path, 'marker_based_metrics.txt'];
%markers_target_filename = [results_path, marker_base_metics_forlder, experiment_name, '.txt'];
%copyfile(markers_source_filename, markers_target_filename);

%% Online metrics
online_metrics_folder = 'online\';
online_source_filename = [sequence_path, 'online_continuous_metrics.txt'];
online_target_filename = [results_path, online_metrics_folder, experiment_name, '.txt'];
copyfile(online_source_filename, online_target_filename);