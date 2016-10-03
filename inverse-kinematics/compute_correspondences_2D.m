function [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points)

centers = cell(length(segments), 1);
for i = 1:length(segments)
    centers{i} = segments{i}.global(1:2, 4);
end
[model_indices, model_points, cell_block_indices] = compute_skeleton_projections(data_points, centers, blocks);

segment_indices = zeros(length(model_indices), 1);
for i = 1:length(model_indices)
    for j = length(blocks):-1:1
        if all(ismember(model_indices{i}, blocks{j}))
            segment_indices(i) = j;
        end
    end
end