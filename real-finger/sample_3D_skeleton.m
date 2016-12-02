function [data_points, segment_indices] = sample_3D_skeleton(segments, num_samples)

B = 3;
data_points = cell(B * num_samples, 1);
segment_indices = zeros(B * num_samples, 1);

k = 1;
for i = 1:length(segments)
    if segments{i}.parent_id == 0, continue; end
    a = segments{i}.global(1:3, 4);
    b = segments{segments{i}.parent_id}.global(1:3, 4);
    x = linspace(a(1), b(1), num_samples);
    y = linspace(a(2), b(2), num_samples);
    z = linspace(a(3), b(3), num_samples);
    for j = 1:num_samples
        data_points{k} = [x(j); y(j); z(j)];   
        segment_indices(k) = i - 1;
        k = k + 1;
    end    
end


