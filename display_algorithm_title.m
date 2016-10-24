function [] = display_algorithm_title(w2, settings, problem_type)

set(gca, 'fontSize', 12); set(gca,'fontname','Cambria');
algorithm_name = '';
if settings.quadratic_one == true, algorithm_name = 'quadratic-one'; end
if settings.quadratic_two == true, algorithm_name = 'quadratic-two'; end
if settings.kalman_like == true, algorithm_name = 'kalman-like'; end
if settings.independent == true, algorithm_name = 'independent'; end

if settings.batch == true && settings.batch_online == true, algorithm_name = 'batch-online'; end
if settings.batch == true && settings.batch_independent == true, algorithm_name = 'batch-independent'; end
if settings.batch == true && settings.batch_online_robust == true, algorithm_name = 'batch-online-robust'; end

title_string = ['\color[rgb]{0.9 0.4 0.3}', algorithm_name, '\color[rgb]{0.25 0.25 0.25}'];

if settings.batch, title_string = [title_string, '     batch = ', '\color[rgb]{0 0.6 0.3}', num2str(settings.batch_size), '\color[rgb]{0.25 0.25 0.25}']; end
if settings.batch && settings.batch_online_robust, title_string = [title_string, '     \tau = ', '\color[rgb]{0 0.6 0.3}', num2str(settings.batch_online_robust_tau), '\color[rgb]{0.25 0.25 0.25}']; end

title_string = [title_string, '     \omega_2 = ', '\color[rgb]{0 0.6 0.3}', num2str(w2), '\color[rgb]{0.25 0.25 0.25}'];
title_string = [title_string, '     \sigma_{data} = ', '\color[rgb]{0 0.6 0.3}', num2str(settings.measurement_noise_std), '\color[rgb]{0.25 0.25 0.25}'];

if strcmp(problem_type, 'sticks_finger')
    title_string = [title_string, '     \sigma_{\beta} = ', '\color[rgb]{0 0.6 0.3}', num2str(settings.beta_noise_std), '\color[rgb]{0.25 0.25 0.25}'];
    title_string = [title_string, '     \sigma_{\theta} = ', '\color[rgb]{0 0.6 0.3}', num2str(settings.theta_noise_std), '\color[rgb]{0.25 0.25 0.25}'];
end

title(title_string, 'interpreter','tex', 'FontWeight','Normal', 'fontsize', 16);

