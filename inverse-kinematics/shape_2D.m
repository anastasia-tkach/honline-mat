function [segments] = shape_2D(segments, beta)

for i = 2:length(segments)
    parent_segment_id = uint16(segments{i}.parent_id);
    segments{i}.local(2, 4) = beta(parent_segment_id);
end

