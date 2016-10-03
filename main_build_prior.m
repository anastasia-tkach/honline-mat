clear; clc; close all;
figure_borders = [0.05 0.08 0.93 0.90];
rng(568);

%% Parameters
B = 3;
T = 3;
beta_noise_std = 0.1;
theta_noise_std = 0.15;
measurement_noise_std = 0.1;

num_experiments = 70;
num_steps = 20;
betas = zeros(num_steps, num_steps, num_experiments, B);
thetas = zeros(num_steps, num_steps, num_experiments, T);
betas_init = zeros(num_steps, num_steps, B);
thetas_init = zeros(num_steps, num_steps, T);

%% Experiments
steps = linspace(-pi/2, pi/2, num_steps);
for t2 = 1:num_steps
    disp(t2);
    for t3 = 1:num_steps
        disp(t3);
        beta_init = [3; 3; 3];
        theta_init = [0; steps(t2); steps(t3)];
        betas_init(t2, t3, :) = beta_init;
        thetas_init(t2, t3, :) = theta_init;
        
        for  i = 1:num_experiments
            
            beta_true = beta_init + beta_noise_std * randn(B, 1);
            theta_true = theta_init + theta_noise_std * randn(T, 1);
            [beta, theta] = fit_one_pose(beta_init, beta_true, theta_init, theta_true, measurement_noise_std, false);
            
            betas(t2, t3, i, :) = beta - beta_true;
            thetas(t2, t3, i, :) = theta - theta_true;
        end
    end
end

%% Display
writing_video = false;
if (writing_video)
    video_writer = VideoWriter('C:\Users\t-antka\Desktop\newfile.avi');
    video_writer.FrameRate = 1; video_writer.Quality = 100; open(video_writer);
end

betas_prior_std = zeros(num_steps, num_steps, B);
for t2 = 1:num_steps
    f = figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.33, 0.8]); hold on;
    set(gca,'position', figure_borders, 'units','normalized'); set(gcf,'color', [1.0, 0.97, 0.96]);

    h = subplot(B, 1, 1); hold on; axis equal; axis off;
    p = get(h, 'pos');
    set(h, 'pos', [0.05, p(2) - 0.06, 0.9, 0.28]);
    for t3 = 1:num_steps
        theta = [thetas_init(t2, t3, 1); thetas_init(t2, t3, 2); thetas_init(t2, t3, 3)];
        display_posed_model(beta_init, theta, [5; 5; 5; 5], [1.0, 0.6, 0.5]);
    end

    for b = 1:B - 1
        h = subplot(B, 1, b + 1); hold on;
        p = get(h, 'pos');
        set(h, 'pos', [0.05, p(2) - 0.06, 0.9, 0.28]);
        %set(h, 'pos', [0.05, p(2) - 0.06, 0.9, 0.42]);
        set(gca,'XTick',[]);
        myline([thetas_init(t2, 1, 3); 0], [thetas_init(t2, num_steps, 3); 0], [1.0, 0.6, 0.5]);
        
        for t3 = 1:num_steps
            theta = thetas_init(t2, t3, 3);
            betas_prior_std(t2, t3, b) = std(betas(t2, t3, :, b));
            myline([theta; betas_prior_std(t2, t3, b)], [theta; -betas_prior_std(t2, t3, b)], [0.65, 0.9, 0.9], 13);
            scatter(theta * ones(num_experiments, 1) + 0.005 * randn(num_experiments, 1),...
                betas(t2, t3, :, b), 25, [0.3, 0.6, 0.8], 'filled');
        end
        ylim([-3, 3]);
        xlim([steps(1)- 0.01, steps(num_steps)+ 0.01]);
    end
    if (writing_video)
        frame = getframe(f);
        writeVideo(video_writer, frame.cdata);
    end
end

% save betas_prior_std betas_prior_std;
% save thetas_init thetas_init;
if exist('video_writer', 'var'), video_writer.close(); end
