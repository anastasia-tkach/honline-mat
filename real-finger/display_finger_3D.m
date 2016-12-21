function [] = display_finger_3D(segments, radii_array, blocks, data_points, model_points)

segment_colors = {[], [255, 173, 153]/255, [179, 220, 160]/255, [255, 173, 153]/255};

%% Display skeleton
%%{
mypoints({segments{1}.global(1:2, 4)}, [255, 173, 153]/255, 35);

for i = 1:length(segments)
    if segments{i}.parent_id == 0, continue; end
    a = segments{i}.global(1:3, 4);
    b = segments{segments{i}.parent_id}.global(1:3, 4);
    myline(a, b, segment_colors{i}, 6);
end
mypoints({a}, [255, 173, 153]/255, 35);
%%}
mylines(data_points, model_points, [0.85, 0.85, 0.85]);
mypoints(model_points, 'b', 10);
mypoints(data_points, 'r', 10);

%% Display surface
centers = cell(length(segments), 1);
radii = cell(length(segments), 1);
for i = 1:length(segments)
    centers{i} = segments{i}.global(1:3, 4);
    radii{i} = radii_array(i);
end
display_result(centers, [], [], blocks(3), radii, false, 0.5, 'none');
view([180, -90]); camlight;

%ylim([-2, 13]); xlim([-7, 7]);
