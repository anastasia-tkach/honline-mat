function [settings, history] = set_batch_size(settings)

B = 3; T = 3;

active_algorithms_count = 0;

if settings.quadratic_one == true
    settings.batch_size = 2;
    active_algorithms_count = active_algorithms_count + 1;
end
if settings.quadratic_two == true
    settings.batch_size = 2;
    active_algorithms_count = active_algorithms_count + 1;
end
if settings.kalman_like == true
    settings.batch_size = 1;
    active_algorithms_count = active_algorithms_count + 1;
end
if settings.kalman_two == true
    settings.batch_size = 2;
    active_algorithms_count = active_algorithms_count + 1;
end
if settings.batch == true
    active_algorithms_count = active_algorithms_count + 1;
end
if settings.independent == true
    settings.batch_size = 1;
    active_algorithms_count = active_algorithms_count + 1;
end
if settings.balman == true
	active_algorithms_count = active_algorithms_count + 1;
end

history.x_batch = zeros(settings.num_frames, (B + T) * settings.batch_size);
history.h_batch = zeros(settings.num_frames, (B + T) * settings.batch_size);

if settings.display_covariance
    history.covariance = zeros(settings.num_frames, B, B);
end

%if settings.balman_kalman_prior || ~settings.balman_solve_all
    history.hessian_independent = zeros(settings.num_frames, B, B);
    history.mu_independent = zeros(settings.num_frames, B);
%end

%% Check if there is only one algorithm set to true

if active_algorithms_count == 0
    error('NO ACTIVE ALGORITHM');
end
if active_algorithms_count > 1
    error('TWO ACTIVE ALGORITHMS');
end

%% Check if only one type of batch is set to true

if settings.batch == true
    if ~settings.batch_independent && ~settings.batch_online && ~settings.batch_online_robust;
        error('BATCH TYPE IS NOT SET');
    end
    if (settings.batch_independent && settings.batch_online) || ...
            (settings.batch_independent && settings.batch_online_robust) || ...
            (settings.batch_online && settings.batch_online_robust)
        error('TWO BATCH TYPES');
    end
end

%% Check if only one type of quadratic is set to true

if settings.quadratic_two == true
    if ~settings.quadratic_two_maximization && ~settings.quadratic_two_marginalization
        error('QUADRATIC TYPE IS NOT SET');
    end
    if settings.quadratic_two_maximization && settings.quadratic_two_marginalization
        error('TWO QUADRATIC TYPES');
    end
end

