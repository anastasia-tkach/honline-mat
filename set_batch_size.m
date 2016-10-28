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
if settings.batch == true
    active_algorithms_count = active_algorithms_count + 1;
end
if settings.independent == true
    settings.batch_size = 1;
    active_algorithms_count = active_algorithms_count + 1;
end

history.x_batch = zeros(settings.num_frames, (B + T) * settings.batch_size);
history.h_batch = zeros(settings.num_frames, (B + T) * settings.batch_size);

if settings.display_covariance
    history.covariance = zeros(settings.num_frames, B, B);
end

%% Check if there is only one algorithm set to true

if active_algorithms_count == 0
    error('NO ACTIVE ALGORITHM');
end
if active_algorithms_count > 1
    error('TWO ACTIVE ALGORITHMS');
end