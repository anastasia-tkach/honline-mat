function [model_points, segment_indices] = compute_correspondences_3D(segments, radii_array, blocks, data_points)

centers = cell(length(segments), 1);
radii = cell(length(segments), 1);
for i = 1:length(segments)
    centers{i} = segments{i}.global(1:3, 4);
    radii{i} = radii_array(i);
end

[model_indices, model_points, block_indices, axis_projections, is_best_projection] = compute_projections_front(data_points, centers, blocks, radii, [0; 0; 1]); 

segment_indices = zeros(length(model_indices), 1);
for i = 1:length(model_indices)
    for j = length(blocks):-1:1
        if all(ismember(model_indices{i}, blocks{j}))
            segment_indices(i) = j;
        end
    end
end