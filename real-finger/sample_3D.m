function [data_points] = sample_3D(segments, radii_array, blocks)

D = 3; 
settings.fov = 15;
downscaling_factor = 80;
settings.H = round(460/downscaling_factor);
settings.W = round(640/downscaling_factor);
settings.D = D;
settings.sparse_data = false;
settings.RAND_MAX = 32767;
settings.side = 'front';
settings.view_axis = 'Z';

%% get centers
centers = cell(length(segments), 1);
radii = cell(length(segments), 1);
for i = 1:length(segments)
    centers{i} = segments{i}.global(1:3, 4);
    radii{i} = radii_array(i);
end

%% Render the data
data_bounding_box = compute_model_bounding_box(centers, radii);

[raytracing_matrix, ~, camera_center] = get_raytracing_matrix(centers, radii, data_bounding_box, settings.view_axis, settings, settings.side);
[rendered_model, I] = render_tracking_model(centers, blocks, radii, raytracing_matrix, camera_center, settings);

[I, J] = find((rendered_model(:, :, 3) > - settings.RAND_MAX));

data_points  = cell(length(I), 1);
for k = 1:length(I),
    data_points{k} = squeeze(rendered_model(I(k), J(k), :));
end
