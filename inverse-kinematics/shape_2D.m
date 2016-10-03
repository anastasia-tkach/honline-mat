function [segments] = shape_2D(segments, beta)

for i = 2:length(segments)
    segments{i}.local(2, 4) = beta(segments{i}.parent_id);
end

