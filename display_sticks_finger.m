function [] = display_sticks_finger(segments, data_points, model_points)

dark_red = [188, 58, 117]/255;
grey = [110, 130, 130]/255;
light_red = [238, 198, 199]/255;
dark_green = [61, 131, 119]/255;
light_green = [144, 194, 171]/230;
orange = [255, 173, 153]/255;

%% Figure 2
% segment_colors = {[], light_red, light_green, light_green};
% mypoints({segments{1}.global(1:2, 4)}, light_red, 600);

%% Figure 3
segment_colors = {[], light_green, light_red, light_green};
mypoints({segments{1}.global(1:2, 4)}, light_green, 600);


mylines(data_points, model_points, [0.85, 0.85, 0.85]);
for i = 1:length(segments)
    if segments{i}.parent_id == 0, continue; end
    a = segments{i}.global(1:2, 4);
    b = segments{segments{i}.parent_id}.global(1:2, 4);
    myline(a, b, segment_colors{i}, 25);
    mypoints({a}, segment_colors{i}, 600);
    mypoints({b}, segment_colors{i}, 600);
end

mypoints(data_points, dark_green, 70, 'o');
ylim([-2, 13]); xlim([-7, 7]);
