function [ellipse_points, eigenval, eigenvec, ok] = get_covarince_elipse(sigma, chisquare_val)

ok = true;

if any(isinf(sigma))
    ok = false;
    ellipse_points = [];
    eigenval = [];
    eigenvec = [];
    return;
end

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
a = chisquare_val*sqrt(largest_eigenval);
b = chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates
ellipse_x_r = a*cos(linspace(0,2*pi));
ellipse_y_r = b*sin(linspace(0,2*pi));

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
ellipse_points = [ellipse_x_r;ellipse_y_r]' * R;