function [F, J, H] = sticks_finger_fg_data(x, segments0, joints, data_points, settings, jacobian_type)

%disp(x);
global video_writer

blocks = {[1, 2], [2, 3], [3, 4]};
B = 3; T = 3;

%% initialize

beta = x(1:B);
theta = x(B + 1:B + T);
[segments] = shape_2D(segments0, beta);
[segments] = pose_2D(segments, joints, theta);


%% data-model correspondences
model_projections = {};
data_samples = {};
model_samples = {};
real_data_points = data_points;

if settings.data_model_energy
    [segment_indices, model_points] = compute_data_correspondences_cpp_wrapper(segments, blocks, data_points);
    model_projections = model_points;
    
    %% single point
    %{
    %data_points = data_points(3);
    %model_points = model_points(3);
    %real_data_points = real_data_points(3);
    %model_projections = model_projections(3);
    %segment_indices = segment_indices(3);
    %}
end


%% model-data correspondences
if (settings.model_data_energy && settings.silhouette_energy)
    error('both model_data_energy and silhouette_energy are active');
end

if settings.model_data_energy || settings.silhouette_energy
    [model_samples, sample_segment_indices] = sample_2D(segments, settings.num_samples);
    
    %% single point
    %{
    data_points = data_points(3);
    model_samples = model_samples(3);
    sample_segment_indices = sample_segment_indices(3);
    real_data_points = real_data_points(3);
    %}
    
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
            if min_distance <= data_radius
                model_samples(k) = [];
                data_samples(k) = [];
                sample_segment_indices(k) = [];
            else
                data_samples{k} = model_samples{k} + (min_distance - data_radius) / min_distance * (data_samples{k} - model_samples{k});
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

if settings.display_iterations
    [segments0, joints] = segments_and_joints_2D();
    [segments0] = shape_2D(segments0, beta);
    [segments] = pose_2D(segments0, joints, theta);
    clf; hold on; axis off; axis equal; set(gcf,'color','w');
    if settings.silhouette_energy
        for i = 1:length(data_samples)
            draw_circle(data_samples{i}, data_radius, [0.85, 0.85, 0.85], 2);
        end
    end
    display_sticks_finger(segments, {}, model_points);
    
    d2m_color = [74, 154, 99]/255;
    m2d_color = [255, 83, 40]/255;
    
    mypoints(real_data_points, d2m_color, 15);
    mypoints(model_projections, d2m_color, 15);
    mylines(real_data_points, model_projections, d2m_color);
    
    mypoints(model_samples, m2d_color, 15);
    mypoints(data_samples, m2d_color, 15);
    mylines(model_samples, data_samples, m2d_color);
    
    drawnow; pause(0.0); %waitforbuttonpress;
    
    if settings.write_video
        f = getframe();
        writeVideo(video_writer, f.cdata);
    end
end

%% Compute Jacobians
[F, J, H] = jacobian_shape_pose_cpp_wrapper(beta, theta, segments, joints, model_points, data_points, segment_indices, jacobian_type);

%% Put weights
if settings.data_model_energy && (settings.model_data_energy || settings.silhouette_energy)
    if ~isempty(data_samples)
        L = size(data_samples, 1);
        F(end - L + 1:end) = sqrt(settings.w1) * F(end - L + 1:end);
        J(end - L + 1:end, :) = sqrt(settings.w1) * J(end - L + 1:end, :);
    end
end


