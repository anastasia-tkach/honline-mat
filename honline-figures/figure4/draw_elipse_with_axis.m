function [] = draw_elipse_with_axis(mu, sigma, dark_color, light_color)

chisquare_val = 2.4477;

[ellipse_points, eigenval, eigenvec] = get_covarince_elipse(sigma, 2.4477);
start1 = mu - chisquare_val * sqrt(eigenval(1)) * eigenvec(:, 1); end1 = mu + chisquare_val * sqrt(eigenval(1)) * eigenvec(:, 1);
start2 = mu - chisquare_val * sqrt(eigenval(2)) * eigenvec(:, 2); end2 = mu + chisquare_val * sqrt(eigenval(2)) * eigenvec(:, 2);
plot([start1(1), end1(1)], [start1(2), end1(2)], 'color', light_color, 'linestyle', '-', 'lineWidth', 2);
plot([start2(1), end2(1)], [start2(2), end2(2)], 'color', light_color, 'linestyle', '-', 'lineWidth', 2);
plot(ellipse_points(:,1) + mu(1), ellipse_points(:,2) + mu(2), '-', 'lineWidth', 3, 'color', dark_color);
mypoint(mu, dark_color, 50);



