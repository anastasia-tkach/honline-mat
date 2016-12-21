%%{
close all;
clear;
%clc;
rng(10);
% load random_seed;
% rng(random_seed);

Rx = @(alpha) [1, 0, 0; 0, cos(alpha), -sin(alpha); 0, sin(alpha), cos(alpha)];
Ry = @(alpha) [cos(alpha), 0, sin(alpha); 0, 1, 0; -sin(alpha), 0, cos(alpha)];
Rz = @(alpha)[cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1];

beta = 5 * rand(1, 1);
t = randn(3, 1);
theta = -pi + 2 * pi * rand(3, 1);
centers_initial{1} = [0; 0; 0];
centers_initial{2} = beta * [0; 1; 0];

while(true)
    c1 = centers_initial{1}; c2 = centers_initial{2};
    x1 = 2 * rand(1, 1); x2 = 0.2 * rand(1, 1); 
    r1 = max(x1, x2); r2  = min(x1, x2);
    if norm(c1 - c2) > r1 + r2
        break;
    end
end
radii = {r1; r2};
blocks = {[1, 2]};


R = @(theta) Rz(theta(1)) * Ry(theta(2)) * Rz(theta(3));

T = [[R(theta), t];[0, 0, 0, 1]];


centers = cell(length(centers_initial), 1);
for i = 1:2
    centers{i} = transform(centers_initial{i}, T);    
end

while(true)
    p = randn(3, 1);
    [index, q, s, is_inside] = projection_convsegment(p, centers{1}, centers{2}, radii{1}, radii{2}, 1, 2);
    if length(index) == 2
        break;
    end
end

frame{1} = [1; 0; 0];
frame{2} = [0; 1; 0];
frame{3} = [0; 0; 1];
for f = 1:3
    rotated_frame{f} = T(1:3, 1:3) * frame{f};
end
%%}
%% compute parameterization
[t] = project_point_on_plane(q, s, rotated_frame{2});
delta_minus = get_agnle_3D(rotated_frame{2}, q - s, []);
delta = pi/2 - delta_minus;
alpha = get_agnle_3D(rotated_frame{1}, t - s, rotated_frame{2});

k_beta = norm(centers{1} - s) / beta;

[rotated_by_alpha] = rotate_around_axis(rotated_frame{2}, rotated_frame{1}, alpha);
[rotated_by_delta] = rotate_around_axis(cross(rotated_by_alpha, rotated_frame{2}), rotated_by_alpha, delta);

[q_, dq] = jacobian_convsegment(q, centers{1}, centers{2}, radii{1}, radii{2}, {'c2'});
dq_analyt = dq.dc2 * rotated_frame{2};

[q_beta, dq_beta] = jacobian_convsegment_beta(q, centers{1}, beta, r1, r2, rotated_frame{2}, {'beta'});

dq_analyt = dq_beta.beta;


%% display

display_result(centers, [], [], blocks, radii, false, 0.5, 'big'); hold on;
for f = 1:2, myline(s, s + rotated_frame{f}, 'k'); end

myline(centers{1}, centers{2}, 'b', 2);
mypoint(q, 'b', 30);
mypoint(s, 'k', 30);

myline(q, t, 'b');
myline(s, t, 'c');

myline(s, s + rotated_by_alpha, [0.7, 0.7, 0.7]);
r = get_convolution_radius_at_points(centers, radii, [1, 2], [], s);
myline(s, s + r * rotated_by_delta, 'm');

%% change parameter 
s0 = s;
offset = norm(centers{1} - s0);

e = 10e-10;
beta = beta + e;
centers{2} = centers{1} + beta * rotated_frame{2};

[index, q_new, s_new, is_inside] = projection_convsegment(q, centers{1}, centers{2}, radii{1}, radii{2}, 1, 2);
delta_minus = get_agnle_3D(rotated_frame{2}, q_new - s_new, []);
delta = pi/2 - delta_minus;

[rotated_by_alpha] = rotate_around_axis(rotated_frame{2}, rotated_frame{1}, alpha);
[rotated_by_delta] = rotate_around_axis(cross(rotated_by_alpha, rotated_frame{2}), rotated_by_alpha, delta);

s = centers{1} + k_beta * beta * rotated_frame{2};

r = get_convolution_radius_at_points(centers, radii, [1, 2], [], s_new);
q_plus = s_new + r * rotated_by_delta;

dq_num = (q_plus - q) / e;

%% display
display_result(centers, [], [], blocks, radii, false, 0.5, 'none'); hold on;

myline(s, q_plus, 'm');
mypoint(q_plus, 'm', 30);
myline(q, q + dq_analyt / norm(dq_analyt), 'g');
myline(q, q + dq_num / norm(dq_num), 'y');
mypoint(s_new, 'c', 30);
%% numerical gradient
disp([dq_analyt, dq_num]);






