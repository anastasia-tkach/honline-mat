function [F, J] = jacobian_shape_2D(segments, model_points, data_points, segment_indices)
num_model_points = size(model_points, 1);
J = zeros(num_model_points, length(segments) - 1);
F = zeros(num_model_points, 1);
u = [0; 1; 0];

%% Build the Jacobian matrix
for k = 1:num_model_points
    d = [data_points{k}; 0];
    m = [model_points{k}; 0];
    
    if isempty(m) || norm(m - d) == 0, continue; end
    n = (d - m) / norm(m - d);
    
    j = zeros(3, length(segments) - 1);
    
    segment = segments{segment_indices(k)};
    for l = 1:length(segment.shape_chain)
        segment_id = segment.shape_chain(l);
        p = segments{segment_id}.global(1:3, 4);
        T = segments{segment_id}.global;
        v = T * [u; 1]; v = v(1:3) / v(4);
        v = v - p;
        
        j(:, segment_id) = v;
    end
    
    %% accumulate sides
    J(k, :) = n' * j;
    F(k) = n' * (d - m);
end
