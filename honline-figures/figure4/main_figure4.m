clc; clear;
n = 2;

x1 = randn(n, 1);
x2 = randn(n, 1);
j1 = randn(n, n);
j2 = randn(n, n);
h1 = j1' * j1;
h2 = j2' * j2;

h1 = 0.5 * [0.1748   -0.1091
   -0.1091    0.15];

h2 = 3 * [0.1748   0.1091
   0.1091    0.075];
    
x1 = [-12.1520; -2.6792];
x2 = [3.0050; -1];

sigma1 = inv(h1);
sigma2 = inv(h2);

x3 = inv(h1 + h2) * (h1 * x1 + h2 * x2);
sigma3 = inv(h1 + h2);

%% display
dark_red = [188, 58, 117]/255;
light_red = [230, 168, 169]/255;
dark_green = [61, 131, 119]/255;
light_green = [144, 194, 171]/240;
light_grey = [195, 195, 195,]/255;
dark_grey = [150, 150, 150]/255;

color2 = dark_green;
color3 = dark_red;

figure; axis off; axis equal; hold on; 
draw_elipse_with_axis(x1, sigma1, dark_grey, light_grey);
draw_elipse_with_axis(x2, sigma2, dark_green, light_green);
draw_elipse_with_axis(x3, sigma3, dark_red, light_red);

ellipse_points = get_covarince_elipse(sigma1, 2.4477);
plot(ellipse_points(:,1) + x1(1), ellipse_points(:,2) + x1(2), '-', 'lineWidth', 3, 'color', dark_grey);

ellipse_points = get_covarince_elipse(sigma2, 2.4477);
plot(ellipse_points(:,1) + x2(1), ellipse_points(:,2) + x2(2), '-', 'lineWidth', 3, 'color', dark_green);

ellipse_points = get_covarince_elipse(sigma3, 2.4477);
plot(ellipse_points(:,1) + x3(1), ellipse_points(:,2) + x3(2), '-', 'lineWidth', 3, 'color', dark_red);
