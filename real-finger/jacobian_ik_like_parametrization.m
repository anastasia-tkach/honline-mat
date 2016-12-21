%%{
close all;
clear;
% clc;
%rng(10);
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

parameter = 'beta';

%% parametrization

k_beta = norm(s - centers{1}) / beta;
offset = inv(R(theta)) * (q - s);

%% analytical jacobian
dq_analyt = k_beta * rotated_frame{2};

%% function-plus
e = 10e-10;

beta_plus = beta + e;
s_plus = centers{1} + k_beta * beta_plus * rotated_frame{2};
q_plus = s_plus + R(theta) * offset;

%% function-minus

beta_minus = beta - e;
s_minus = centers{1} + k_beta * beta_minus * rotated_frame{2};
q_minus = s_minus + R(theta) * offset;

%% numerical jacobian
dq_num = (q_plus - q_minus) / 2 / e;

%% display
%{
display_result(centers, [], [], blocks, radii, false, 0.5, 'big'); hold on; 
mypoint(q, 'k', 30);
myline(q, q + dq_analyt, 'b');
myline(q, q + dq_num, 'r');

display_result(centers, [], [], blocks, radii, false, 0.5, 'none'); hold on; 
mypoint(q_plus, 'r', 30);
%}

disp([dq_analyt'; dq_num']);


