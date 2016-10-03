function [data_points] = sample_2D(segments, num_samples)

data_points = {};

for i = 1:length(segments)
    if segments{i}.parent_id == 0, continue; end
    a = segments{i}.global(1:2, 4);
    b = segments{segments{i}.parent_id}.global(1:2, 4);
    x = linspace(a(1), b(1), num_samples);
    y = linspace(a(2), b(2), num_samples);
    for j = 1:num_samples
        data_points{end + 1} = [x(j); y(j)];  %#ok<AGROW>        
    end    
end


