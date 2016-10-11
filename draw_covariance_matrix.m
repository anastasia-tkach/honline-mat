function [] = draw_covariance_matrix(mu, sigma)

[eigenvec, eigenval] = eig(sigma);
eigenval = diag(eigenval);
largest_eigenval = max(eigenval);
largest_eigenvec = eigenvec(:, eigenval == largest_eigenval);
smallest_eigenval = min(eigenval);
smallest_eigenvec = eigenvec(:, eigenval == smallest_eigenval);

% Calculate the angle between the x-axis and the largest eigenvector
phi = atan2(largest_eigenvec(2), largest_eigenvec(1));
if(phi < 0), phi = phi + 2*pi; end

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
a = chisquare_val*sqrt(largest_eigenval);
b = chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r = a*cos(linspace(0,2*pi));
ellipse_y_r = b*sin(linspace(0,2*pi));

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

%% Draw the error ellipse
figure; axis equal; hold on;
plot(r_ellipse(:,1) + mu(1), r_ellipse(:,2) + mu(2), '-', 'lineWidth', 2, 'color', [0.1, 0.5, 0.7]);
myline(mu, mu + chisquare_val * sqrt(eigenval(1)) * eigenvec(:, 1), [0.7, 0, 0.3]);
myline(mu, mu + chisquare_val * sqrt(eigenval(2)) * eigenvec(:, 2), [1, 0.4, 0]);
mypoint(mu, [0, 0.5, 0.7], 50);
legend({'elipse', 'x_{n-1}', 'x_n'})

%% Explore

A = sigma(1, 1);
B = sigma(1, 2);
D = sigma(2, 2);
disp(D);

myline(mu - chisquare_val * sqrt(A) * [1; 0] - chisquare_val * sqrt(D) * [0; 1],...
    mu - chisquare_val * sqrt(A) * [1; 0] + chisquare_val * sqrt(D) * [0; 1], [0.7, 0.7, 0.7]);
myline(mu - chisquare_val * sqrt(A) * [1; 0] - chisquare_val * sqrt(D) * [0; 1],...
    mu + chisquare_val * sqrt(A) * [1; 0] - chisquare_val * sqrt(D) * [0; 1], [0.7, 0.7, 0.7]);
myline(mu + chisquare_val * sqrt(A) * [1; 0] + chisquare_val * sqrt(D) * [0; 1],...
    mu - chisquare_val * sqrt(A) * [1; 0] + chisquare_val * sqrt(D) * [0; 1], [0.7, 0.7, 0.7]);
myline(mu + chisquare_val * sqrt(A) * [1; 0] + chisquare_val * sqrt(D) * [0; 1],...
    mu + chisquare_val * sqrt(A) * [1; 0] - chisquare_val * sqrt(D) * [0; 1], [0.7, 0.7, 0.7]);


xlim([mu(1) - chisquare_val * sqrt(A) - 0.5, mu(1) + chisquare_val * sqrt(A) + 0.5]);
ylim([mu(2) - chisquare_val * sqrt(D) - 0.5, mu(2) + chisquare_val * sqrt(D) + 0.5]);
