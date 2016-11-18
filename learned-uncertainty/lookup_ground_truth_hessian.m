function [h, theta_closest] = lookup_ground_truth_hessian(theta_lookup, thetas, true_hessians)

thetas2_index = thetas(:, 1, 2);
thetas3_index = thetas(1, :, 3);
[~, lookup_index_theta2] = min(abs(thetas2_index - theta_lookup(2)));
[~, lookup_index_theta3] = min(abs(thetas3_index - theta_lookup(3)));
h = squeeze(true_hessians(lookup_index_theta2, lookup_index_theta3, :, :));

theta_closest = squeeze(thetas(lookup_index_theta2, lookup_index_theta3, :));