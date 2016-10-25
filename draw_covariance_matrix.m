function [] = draw_covariance_matrix(mu, sigma, frame_certainty)

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
mu = [0; 0];

uncertain_background_color = [1; 0.98; 0.95];
certain_background_color = [0.96; 1; 0.93];

figure; axis equal; hold on;
if frame_certainty
    set(gca,'color', certain_background_color);
else
    set(gca,'color', uncertain_background_color);
end
plot(r_ellipse(:,1) + mu(1), r_ellipse(:,2) + mu(2), '-', 'lineWidth', 2, 'color', [0.1, 0.5, 0.7]);
myline(mu, mu + chisquare_val * sqrt(eigenval(1)) * eigenvec(:, 1), [0.7, 0, 0.3]);
myline(mu, mu + chisquare_val * sqrt(eigenval(2)) * eigenvec(:, 2), [1, 0.4, 0]);


%% Explore

A = sigma(1, 1);
B = sigma(1, 2);
D = sigma(2, 2);

h = inv(sigma);
%disp(D);

myline(mu - chisquare_val * sqrt(A) * [1; 0] - chisquare_val * sqrt(D) * [0; 1],...
    mu - chisquare_val * sqrt(A) * [1; 0] + chisquare_val * sqrt(D) * [0; 1], [0.1, 0.7, 0.5], 2);

myline( mu - chisquare_val * 1/sqrt(h(1, 1)) * [1; 0] - chisquare_val * 1/sqrt(h(2, 2)) * [0; 1],  ...
    mu - chisquare_val * 1/sqrt(h(1, 1)) * [1; 0] + chisquare_val * 1/sqrt(h(2, 2)) * [0; 1], [0.6, 0.7, 0], 2);

myline(mu - chisquare_val * sqrt(A) * [1; 0] - chisquare_val * sqrt(D) * [0; 1],...
    mu + chisquare_val * sqrt(A) * [1; 0] - chisquare_val * sqrt(D) * [0; 1], [0.7, 0.7, 0.7], 1);
myline(mu + chisquare_val * sqrt(A) * [1; 0] + chisquare_val * sqrt(D) * [0; 1],...
    mu - chisquare_val * sqrt(A) * [1; 0] + chisquare_val * sqrt(D) * [0; 1], [0.7, 0.7, 0.7], 1);
myline(mu + chisquare_val * sqrt(A) * [1; 0] + chisquare_val * sqrt(D) * [0; 1],...
    mu + chisquare_val * sqrt(A) * [1; 0] - chisquare_val * sqrt(D) * [0; 1], [0.1, 0.7, 0.5], 2);


xlim([mu(1) - chisquare_val * sqrt(A) - 0.5, mu(1) + chisquare_val * sqrt(A) + 0.5]);
ylim([mu(2) - chisquare_val * sqrt(D) - 0.5, mu(2) + chisquare_val * sqrt(D) + 0.5]);

%% Display importance


myline( mu - chisquare_val * 1/sqrt(h(1, 1)) * [1; 0] - chisquare_val * 1/sqrt(h(2, 2)) * [0; 1],  ...
    mu + chisquare_val * 1/sqrt(h(1, 1)) * [1; 0] - chisquare_val * 1/sqrt(h(2, 2)) * [0; 1], [0.7, 0.7, 0.7], 1);

myline( mu + chisquare_val * 1/sqrt(h(1, 1)) * [1; 0] + chisquare_val * 1/sqrt(h(2, 2)) * [0; 1],  ...
    mu - chisquare_val * 1/sqrt(h(1, 1)) * [1; 0] + chisquare_val * 1/sqrt(h(2, 2)) * [0; 1], [0.7, 0.7, 0.7], 1);

myline( mu + chisquare_val * 1/sqrt(h(1, 1)) * [1; 0] + chisquare_val * 1/sqrt(h(2, 2)) * [0; 1],  ...
    mu + chisquare_val * 1/sqrt(h(1, 1)) * [1; 0] - chisquare_val * 1/sqrt(h(2, 2)) * [0; 1], [0.6, 0.7, 0], 2);



myline( mu - chisquare_val * 1/sqrt(h(1, 1)) * [1; 0], mu + chisquare_val * 1/sqrt(h(1, 1)) * [1; 0], [0.7, 0.7, 0.7], 1);
myline( mu - chisquare_val * 1/sqrt(h(2, 2)) * [0; 1], mu + chisquare_val * 1/sqrt(h(2, 2)) * [0; 1], [0.7, 0.7, 0.7], 1);

mypoint(mu, [0, 0.5, 0.7], 50);
%lgnd = legend({'3\sigma elipse', 'x_{n-1}', 'x_n', 'marginalization', 'maximization'}); set(lgnd,'color','none', 'edgecolor', 'none');
%xlabel('x_{n-1}');
%ylabel('x_{n}');
xlabel('\beta_1');
ylabel('\beta_2');
set(gca, 'fontSize', 13); set(gca,'fontname','Cambria');




