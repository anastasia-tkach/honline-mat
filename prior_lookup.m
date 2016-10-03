function [beta_std] = prior_lookup(theta, betas_prior_std, thetas_init)
% close all;
% load betas_prior_std;
% load thetas_init;
% T = 3;
% theta_magnitude = pi/2;
% theta = [0; - theta_magnitude + 2 * theta_magnitude * rand(2, 1)];

thetas2_index = thetas_init(:, 1, 2);
thetas3_index = thetas_init(1, :, 3);
[~, lookup_index_theta2] = min(abs(thetas2_index - theta(2)));
[~, lookup_index_theta3] = min(abs(thetas3_index - theta(3)));
beta_std = betas_prior_std(lookup_index_theta2, lookup_index_theta3, :);

return
%% Display
lookup_theta_init = thetas_init(lookup_index_theta2, lookup_index_theta3, :);
beta_init = [3; 3; 3];
figure; axis equal; axis off; hold on;
display_posed_model(beta_init, theta);
xlim([-10, 10]); ylim([-5, 10]);
figure; axis equal; axis off; hold on;
display_posed_model(beta_init, [lookup_theta_init(1); lookup_theta_init(2); lookup_theta_init(3)]);
xlim([-10, 10]); ylim([-5, 10]);

