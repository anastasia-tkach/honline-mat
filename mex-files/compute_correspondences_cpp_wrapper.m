function [segment_indices, model_points] = compute_correspondences_cpp_wrapper(segments, blocks, data_points)

centers = cell(length(segments), 1);
for j = 1:length(segments)
    centers{j} = segments{j}.global(1:2, 4);
end

Centers = zeros(length(centers), 2);
for j = 1:length(centers)
    Centers(j, :) = centers{j}';
end
Blocks = zeros(length(blocks), 2);
for j = 1:length(blocks)
    Blocks(j, :) = blocks{j}';
end
DataPoints = zeros(length(data_points), 2);
for j = 1:length(data_points)
    DataPoints(j, :) = data_points{j}';
end

%[segment_indices, ModelPoints] = compute_correspondences_matlab(Centers, Blocks, DataPoints);

[segment_indices, ModelPoints] = compute_correspondences_cpp([size(Centers, 1), size(Blocks, 1), size(DataPoints, 1)], Centers, Blocks, DataPoints);

model_points = cell(length(data_points), 1);
for j = 1:length(data_points)
    model_points{j} = ModelPoints(j, :)';
end
