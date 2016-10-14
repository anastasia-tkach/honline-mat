clear; close all;

x = -2:0.01:2;

t = 1;
least_squares = x.^2;
german_mccure = @(t) (x/t).^2./(1 + (x/t).^2);

figure; hold on; 
plot(x, least_squares, 'lineWidth', 2);
%plot(x, german_mccure(0.1), 'lineWidth', 2);
%plot(x, german_mccure(0.5), 'lineWidth', 2)
plot(x, german_mccure(0.5), 'lineWidth', 2);
%plot(x, german_mccure(2), 'lineWidth', 2);
ylim([0, 2]);
