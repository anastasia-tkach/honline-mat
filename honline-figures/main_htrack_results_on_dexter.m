data_root = 'E:\OneDrive\EPFL\Code\honline-mat\honline-figures\htrack_results_on_dexter\';

sequence_names = {'adbadd', 'flexex1', 'pinch', 'fingercount', 'tigergrasp','fingerwave', 'random'};

htrack_markers = [];
dexter_markers = [];
num_joints = 24;


for sequence_id = 1:length(sequence_names)
    
    %% Read htrack markers
    data =  load([data_root, 'htrack\', sequence_names{sequence_id}, '.txt']);
    result = zeros(length(data), 5 * 3);
    for i = 1:length(data)
        d = data(i, :);
        d = reshape(d, 3, num_joints)';
        d = d([8 24 20 16 12], :);
        result(i, :) = d(:)';
    end
    htrack_markers = [htrack_markers; result];
    
    %% Load dexter markers
    data =  load([data_root, 'sridhar\', sequence_names{sequence_id}, '.txt']);
    result = zeros(length(data)/6, 6 * 3);
    for i = 1:length(data) / 6        
        d = data(6 * (i - 1) + 1:6 * i , :);
        result(i, :) = d(:)';
    end
    dexter_markers = [dexter_markers; result];
    dexter_markers = dexter_markers(1:size(htrack_markers, 1), :);
end

num_frames = size(htrack_markers, 1);
errors = zeros(num_frames, 5);

% dexter_markers = load('E:\OneDrive\EPFL\Code\honline-mat\honline-figures\htrack_results_on_dexter\sridhar_all.txt');

for i = 1:num_frames
    h = htrack_markers(i, :);
    d = dexter_markers(i, :);
    h = reshape(h, 5, 3);
    d = reshape(d, 6, 3);
    d(:, 1) = -d(:, 1);
    d = d(1:5, :);
    for j = 1:5
        errors(i, j) = norm(h(j, :) - d(j, :));
    end
end

%% Write errors
fileID = fopen('E:\OneDrive\EPFL\Code\honline-mat\honline-figures\htrack_results_on_dexter\errors_all.txt','w');
for i = 1:size(errors, 1)
     for j = 1:size(errors, 2)
        fprintf(fileID,'%8.4f ', errors(i, j));
    end
    fprintf(fileID,'\n');
end
