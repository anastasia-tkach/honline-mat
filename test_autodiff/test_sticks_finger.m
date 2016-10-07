clear; clc;
rng default;
num_samples = 5;
B = 3; T = 3; D = 3 * num_samples;
measurement_noise_std = 0.07;
beta_noise_std = 0.5;
theta_noise_std = 0.15;
blocks = {[1, 2], [2, 3], [3, 4]};
beta = [3; 3; 3];
theta = [pi/3; pi/3; pi/3];
[segments0, joints] = segments_and_joints_2D();

theta_true = theta + theta_noise_std * randn(size(theta));
beta_true = beta + beta_noise_std * randn(size(beta));
[data_segments] = shape_2D(segments0, beta_true);
[data_segments] = pose_2D(data_segments, joints, theta_true);
[data_points] = sample_2D(data_segments, num_samples);
for i = 1:length(data_points)
    data_points{i} = data_points{i} + measurement_noise_std * randn(2, 1);
end

[F_, J_, H_] = sticks_finger_eg_single([beta; theta], segments0, joints, data_points);

%% Change data

%{
segments = segments0;

num_points = length(data_points);
num_segments = length(segments);
num_joints = length(joints);
max_kinematic_chain = 3;
%}
% data
DataPoints = zeros(15, 2);
for j = 1:length(data_points)
    DataPoints(j, :) = data_points{j}';
end
%{
% segments
SegmentsLocal =  zeros(num_segments, 16);
for j = 1:length(segments)
    SegmentsLocal(j, :) = segments{j}.local(:)';
end
SegmentsKinematicChain =  -1 * ones(num_segments, max_kinematic_chain);
for j = 1:length(segments) - 1
    SegmentsKinematicChain(j, 1:length(segments{j}.kinematic_chain)) = segments{j}.kinematic_chain';
end
SegmentsParentId = zeros(num_segments, 1);
for j = 1:num_segments
    SegmentsParentId(j) = segments{j}.parent_id;
end

% joints
JointsSegmentId = zeros(num_joints, 1);
for j = 1:num_joints
    JointsSegmentId(j) = joints{j}.segment_id;
end
JointsAxis = zeros(num_joints, 3);
for j = 1:num_joints
    JointsAxis(j, :) = joints{j}.axis';
end
%}
data = [DataPoints(:)];
    %; SegmentsLocal(:); SegmentsKinematicChain(:); SegmentsParentId(:); JointsSegmentId(:); JointsAxis(:)];

x = [beta; theta];
[F, J, H] = sticks_finger_egf_autodiff(x, data);

au_autodiff_generate(@sticks_finger_egf_autodiff, x, data, 'E:/OneDrive/EPFL/Code/honline-mat/test_autodiff/sticks_finger_egf_cpp.cpp', 'HESSIAN=1');

[V] = au_autodiff_test(x, data, 2);
j_ = V(1:m);
h_ = V(m + 1:end - 1);

J_ = j_';
H_ = zeros(m, m);
k = 0;
for  i = m:-1:1
    H_(1:i, i) = h_(k + 1:k + i);
    H_(i, 1:i-1) =  h_(k + 1:k + i - 1);
    k = k + i;
end

