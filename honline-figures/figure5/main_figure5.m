clear; clc;
%close all;
rng('default')
%% Parameters
settings.num_samples = 10;
B = 3;
T = 3;
settings.measurement_noise_std = 0.001;
settings.beta_bias = [0; 0; 0];
settings.beta_noise_std = 0.0;
settings.theta_noise_std = 0.0;
blocks = {[1, 2], [2, 3], [3, 4]};
[segments0, joints] = segments_and_joints_2D();

settings.data_model_energy = true;
settings.model_data_energy = false;
settings.silhouette_energy = false;
settings.display_converged = false;
settings.write_video = false;
settings.display_iterations = false;
settings.display_jacobian = false;

beta_true = [3; 3; 3];
theta_true = [0;-pi/3; 0]; %-pi/15, -pi/3
beta = beta_true + settings.beta_noise_std * randn;
theta = theta_true + settings.theta_noise_std * randn;
X = [beta; theta];

settings.num_frames = 1;
settings.num_iters = 50;

num_steps = 9;

theta1 = linspace(0.12, pi/2 - 0.12, num_steps);
theta2 = linspace(0.12, pi/2 - 0.12, num_steps);

covariances = cell(num_steps, num_steps);
thetas = cell(num_steps, num_steps);

for t1 = 1:num_steps
    for t2 = 1:num_steps
        %% Initialize
        theta_true = [0; theta1(t1); theta2(t2)];
        [data_segments] = shape_2D(segments0, beta_true);
        [data_segments] = pose_2D(data_segments, joints, theta_true);
        [data_points] = sample_2D(data_segments, settings.num_samples);
        for i = 1:length(data_points)
            data_points{i} = data_points{i} + settings.measurement_noise_std * randn(2, 1);
        end
        
        %% Optimize
        
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_data(X, segments0, joints, data_points, settings, 'cpp'), X, settings.num_iters);
        H = J' * J;
        covariances{t1, t2} = inv(H(1:2, 1:2));
        thetas{t1, t2} = [theta1(t1), theta2(t2)];
        
        %% Display
        %display_jacobian(X, J, settings, N);
        
        if (settings.display_converged)
            figure('units', 'normalized', 'outerposition', [0.25, 0.275, 0.45, 0.7]);
            beta = X(end - (B + T) + 1:end - T);
            beta = [3, 3, 3];
            theta = X(end - T + 1:end);
            data_points = data_points;
            [segments0, joints] = segments_and_joints_2D();
            [segments0] = shape_2D(segments0, beta);
            [segments] = pose_2D(segments0, joints, theta);
            [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
            clf; hold on; axis off; axis equal; set(gcf,'color','w');
            display_sticks_finger(segments, data_points, model_points);
            drawnow; pause(0.05); %waitforbuttonpress
            
        end
    end
end

%% Display energy as a function of beta1
dark_red = [188, 58, 117]/255;
light_red = 0.5 * [217, 154, 143]/255 + 0.5 * [238, 198, 199]/255;
dark_green = [33, 159, 126]/255;
light_green = [144, 194, 171]/230;
orange = [255, 119, 51] / 255;

figure('units', 'normalized', 'outerposition', [0.25, 0.275, 0.30, 0.7]); hold on; axis on; axis equal;
for i = 1:num_steps
    for j = 1:num_steps
        sigma = covariances{i, j};
        mu = thetas{i, j};
        ellipse_points = get_covarince_elipse(sigma, 0.085);
        w1 = norm([i, j]) / norm([num_steps, num_steps]);
        w2 = 1 - w1;
        plot(ellipse_points(:,1) + mu(1), ellipse_points(:,2) + mu(2), '-', 'lineWidth', 3.5, 'color', w1 * light_red + w2 * light_green);
        mypoint(mu, w1 * dark_red + w2 * dark_green, 15);
    end
end
xlim([0, pi/2]);
ylim([0, pi/2]);
box on; set(gca,'linewidth', 1.5); set(gca, 'fontsize', 14, 'fontname', 'Cambria');
set(gca, 'xcolor', [0.5, 0.5, 0.5]);
set(gca, 'ycolor', [0.5, 0.5, 0.5]); 





