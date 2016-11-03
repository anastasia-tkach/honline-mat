function [] = display_jacobian(X, J, settings, N)

if ~settings.display_jacobian
    return
end

B = 3; T = 3;

beta_indices = repmat([ones(B, 1); zeros(B, 1)], min(N, settings.batch_size), 1);

J1 = J(1:min(N, settings.batch_size) * settings.num_samples * 3, :);
JJ1 = J1' * J1;
JJ1 = JJ1(beta_indices == 1, beta_indices == 1);

if (settings.batch || settings.quadratic_one)
    if (N == 2  || (N > 1 && settings.batch_independent))
        J2 = J(min(N, settings.batch_size) * settings.num_samples * 3 + 1:min(N, settings.batch_size) * settings.num_samples * 3 + 3, :);
        JJ2 = J2' * J2;
        JJ2 = JJ2(beta_indices == 1, beta_indices == 1);
        %figure; imagesc(JJ2); axis equal; colorbar;
    end
    if (N > 2 && ~settings.batch_independent)
        J2 = J(min(N, settings.batch_size) * settings.num_samples * 3 + 1:min(N, settings.batch_size) * settings.num_samples * 3 + 6, :);
        JJ2 = J2' * J2;
        JJ2 = JJ2(beta_indices == 1, beta_indices == 1);
        %figure; imagesc(JJ2); axis equal; colorbar;
    end
end
if (settings.uniform_shape_prior)
    J4 = J(min(N, settings.batch_size) * settings.num_samples * 3 + 7:end, :);
    JJ4 = J4' * J4;
    JJ4 = JJ4(beta_indices == 1, beta_indices == 1);
    %figure; imagesc(JJ4); axis equal; colorbar;
end

JJ = J' * J;
JJ = JJ(beta_indices == 1, beta_indices == 1);
figure; imagesc(JJ); axis equal; colorbar; %caxis([cmin cmax]);

mu = X(end - B - T + 1:end - 1 - T);
sigma = inv(JJ);
[ellipse_points] = get_covarince_elipse(sigma(1:2, 1:2), 2.4477);
figure; plot(ellipse_points(:,1) + mu(1), ellipse_points(:,2) + mu(2), '-', 'lineWidth', 2, 'color',  [136, 187, 119]/255); axis equal;