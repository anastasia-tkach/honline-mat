clear; close all;

x = -2:0.01:2;

t = 1;
least_squares = x.^2;
german_mccure = @(t) (x/t).^2./(1 + (x/t).^2);

figure; hold on; 
plot(x, least_squares, 'lineWidth', 2);
plot(x, german_mccure(0.5), 'lineWidth', 2);
plot(x, german_mccure(1.0), 'lineWidth', 2)
plot(x, german_mccure(1.5), 'lineWidth', 2);
plot(x, german_mccure(2.0), 'lineWidth', 2);
ylim([0, 2]);
legend({'least-squares', 't = 0.5', 't = 1.0', 't = 1.5', 't = 2'});
xlabel('x');
ylabel('(x/t)^2/(1 + (x/t))^2');
