function [F, J, H] = sticks_finger_fg_data(x, segments0, joints, data_points, settings)

%disp(x);

blocks = {[1, 2], [2, 3], [3, 4]};
B = 3; T = 3;

%% initialize

beta = x(1:B);
theta = x(B + 1:B + T);
[segments] = shape_2D(segments0, beta);
[segments] = pose_2D(segments, joints, theta);

%% data-model correspondences

if settings.data_model_energy
    [segment_indices, model_points] = compute_data_correspondences_cpp_wrapper(segments, blocks, data_points);
end

%% model-data correspondences
if (settings.model_data_energy && settings.silhouette_energy)
    error('both model_data_energy and silhouette_energy are active');
end

if settings.model_data_energy || settings.silhouette_energy
    [model_samples, sample_segment_indices] = sample_2D(segments, 10);
    
    DataPoints = zeros(length(data_points), 2);
    for j = 1:length(data_points)
        DataPoints(j, :) = data_points{j}';
    end
    ModelSamples = zeros(length(model_samples), 2);
    for j = 1:length(model_samples)
        ModelSamples(j, :) = model_samples{j}';
    end
    D = pdist2(ModelSamples, DataPoints, 'euclidean');
    [~, closest_indices] = min(D, [], 2);
    data_samples = cell(length(model_samples), 1);
    for j = 1:length(closest_indices)
        data_samples{j} = DataPoints(closest_indices(j), :)';
    end
    
    %% Remove the points that are within data silhouette
    if settings.silhouette_energy
        data_radius =  sqrt(settings.measurement_noise_std^2 + 3^2 / settings.num_samples^2);
        k = 1;
        while k <= length(model_samples)
            min_distance = norm(model_samples{k} - data_samples{k});
            if min_distance < data_radius
                model_samples(k) = [];
                data_samples(k) = [];
                sample_segment_indices(k) = [];
            else
                k = k + 1;
            end
        end
        if isempty(sample_segment_indices), sample_segment_indices = sample_segment_indices'; end
    end
end

%% stack together

if settings.data_model_energy && (settings.model_data_energy || settings.silhouette_energy)
    model_points = [model_points; model_samples];
    data_points = [data_points; data_samples];
    segment_indices = [segment_indices; sample_segment_indices];
end

if ~settings.data_model_energy && (settings.model_data_energy || settings.silhouette_energy)
    model_points = model_samples;
    data_points = data_samples;
    segment_indices = sample_segment_indices;
end

%% display
%{
%return
[segments0, joints] = segments_and_joints_2D();
[segments0] = shape_2D(segments0, beta);
[segments] = pose_2D(segments0, joints, theta);
clf; hold on; axis off; axis equal; set(gcf,'color','w');
for i = 1:length(data_points)
    draw_circle(data_points{i}, data_radius, [0.85, 0.85, 0.85], 2);
end
display_sticks_finger(segments, data_points, model_points);
%disp(' ');
waitforbuttonpress
%}

%% Compute Jacobians
[F, J, H] = jacobian_shape_pose_cpp_wrapper(beta, theta, segments, joints, model_points, data_points, segment_indices, 'cpp');


