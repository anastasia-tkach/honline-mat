function [] = display_sticks_finger(segments, data_points, model_points)

segment_colors = {[], [255, 173, 153]/255, [179, 220, 160]/255, [255, 173, 153]/255};

mypoints({segments{1}.global(1:2, 4)}, [255, 173, 153]/255, 35);
mylines(data_points, model_points, [0.85, 0.85, 0.85]);
for i = 1:length(segments)
    if segments{i}.parent_id == 0, continue; end
    a = segments{i}.global(1:2, 4);
    b = segments{segments{i}.parent_id}.global(1:2, 4);
    myline(a, b, segment_colors{i}, 6);
end
mypoints({a}, [255, 173, 153]/255, 35);

mypoints(data_points, [74, 154, 99]/255, 15);
ylim([-2, 13]); xlim([-7, 7]);
drawnow;