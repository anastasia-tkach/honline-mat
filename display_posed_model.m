function [] = display_posed_model(beta, theta, line_width, line_color)

[segments0, joints] = segments_and_joints_2D();

[shaped_segments] = shape_2D(segments0, beta);
[segments] = pose_2D(shaped_segments, joints, theta);

for i = 1:length(segments)
    if segments{i}.parent_id == 0, continue; end
    a = segments{i}.global(1:2, 4);
    b = segments{segments{i}.parent_id}.global(1:2, 4);
    myline(a, b, line_color, line_width(i - 1));
    mypoints({a}, line_color, 30);
    mypoints({b}, line_color, 30);
end