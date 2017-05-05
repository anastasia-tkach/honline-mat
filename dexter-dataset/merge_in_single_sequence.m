
clear;
data_root = 'E:\Data\datasets\dexter1\data\';
sequence_names = {'adbadd', 'flexex1', 'pinch', 'fingercount', 'tigergrasp','fingerwave', 'random'};

Labels_3D = [];
Labels_2D = [];
total_frame_count = 1;
figure;
for sequence_id = 1:length(sequence_names)
    
    %% Read labels
    file_id = fopen([data_root, sequence_names{sequence_id}, '\annotations\Pos3D.txt'],'r');
    labels_3D = fscanf(file_id, '%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n\n');
    labels_3D = reshape(labels_3D, 18, length(labels_3D) / 18)';
    Labels_3D = [Labels_3D; labels_3D];    
    
    file_id = fopen([data_root, sequence_names{sequence_id}, '\annotations\Pos2D.txt'],'r');
    labels_2D = fscanf(file_id, '%f, %f\n%f, %f\n%f, %f\n%f,%f\n%f, %f\n%f,%f\n\n');
    labels_2D = reshape(labels_2D, 12, length(labels_2D) / 12)';
    Labels_2D = [Labels_2D; labels_2D];
    
    %% Read frames
    frame_names = dir([data_root, sequence_names{sequence_id}, '\tof\depth\']);
    frame_names = frame_names(3:end);
    
    for frame_id = 1:length(frame_names)        
        depth = imread([data_root, sequence_names{sequence_id}, '\tof\depth\', frame_names(frame_id).name]);   
        
        %% Show labels
        %%{
        %imshow(depth); hold on;
        %l = Labels_2D(total_frame_count, :);
        %l = reshape(l, 2, 6)';
        %scatter(l(:, 1), l(:, 2), 10, 'y', 'filled');
        %drawnow;
        %%}
        
        %% Write frame
        %{
        filename_prefix = sprintf('%04d', total_frame_count);
        imwrite(depth, [data_root, 'ALL\depth\', filename_prefix, '.png']);
        %}
        total_frame_count = total_frame_count + 1;
    end
    
    Labels_3D = Labels_3D(1:total_frame_count - 1, :);
    Labels_2D = Labels_2D(1:total_frame_count - 1, :);
    disp(sequence_names{sequence_id});
    disp([total_frame_count - 1, size(Labels_3D, 1)]);
    disp('');
end

%% Write 3D labels
%{
fileID = fopen('E:\Data\datasets\dexter1\data\ALL\marker_positions.txt','w');
fprintf(fileID, '%d\n', size(Labels_3D, 1));

for i = 1:size(Labels_3D, 1)
    for j = 1:size(Labels_3D, 2)
        fprintf(fileID,'%8.4f ', Labels_3D(i, j));
    end
    fprintf(fileID,'\n');
end

fclose(fileID);
%}

%% Write 2D labels
%{
Labels_2D = [
    Labels_2D(:, 1:2), Labels_3D(:, 3),...
    Labels_2D(:, 3:4), Labels_3D(:, 6),...
    Labels_2D(:, 4:6), Labels_3D(:, 9),...
    Labels_2D(:, 7:8), Labels_3D(:, 12),...
    Labels_2D(:, 9:10), Labels_3D(:, 15),...
    Labels_2D(:, 11:12), Labels_3D(:, 18)
    ];

fileID = fopen('E:\Data\datasets\dexter1\data\ALL\marker_positions.txt','w');
fprintf(fileID, '%d\n', size(Labels_2D, 1));

for i = 1:size(Labels_2D, 1)
     for j = 1:size(Labels_2D, 2)
        fprintf(fileID,'%8.4f ', Labels_2D(i, j));
    end
    fprintf(fileID,'\n');
end

fclose(fileID);
%}



