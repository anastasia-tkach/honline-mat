clc; clear;
n = 2;

x1 = randn(n, 1);
x2 = randn(n, 1);
j1 = randn(n, n);
j2 = randn(n, n);
h1 = j1' * j1;
h2 = j2' * j2;
sigma1 = inv(h1);
sigma2 = inv(h2);

%% display
color1 = [136, 187, 119]/255;
color2 = [1.0 0.45 0.3];
color3 = [0.7, 0.7, 0.7];
figure; hold on; axis equal;

ellipse_points1 = get_covarince_elipse(sigma1, 1);
plot(ellipse_points1(:,1) + x1(1), ellipse_points1(:,2) + x1(2), '-', 'lineWidth', 2, 'color',  color1);
mypoint(x1, color1, 30);

ellipse_points2 = get_covarince_elipse(sigma2, 1);
plot(ellipse_points2(:,1) + x2(1), ellipse_points2(:,2) + x2(2), '-', 'lineWidth', 2, 'color',  color2);
mypoint(x2, color2, 30);

%% the objective function 
F = @ (x) ... 
    [sqrtm(h1) * (x - x1); ...
    sqrtm(h2) * (x - x2)];
options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt');
x_opt = lsqnonlin(F, x1, [], [], options);

disp(x_opt);

%% closed form
x3_closed_form =  inv(h1 + h2) * (h1 * x1 + h2 * x2);

%% Display
x3 = x_opt;
sigma3 = inv(h1 + h2);
ellipse_points3 = get_covarince_elipse(sigma3, 1);
plot(ellipse_points3(:,1) + x3(1), ellipse_points3(:,2) + x3(2), '-', 'lineWidth', 2, 'color',  color3);
mypoint(x3, color3, 30);
mypoint(x3_closed_form, 'r', 30);


