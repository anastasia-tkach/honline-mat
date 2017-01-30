d = 0:0.01:01;

x = d;

figure; hold on;
plot(d, x, 'linewidth', 2);

y = d - 5;
y(y < 0) = 0;
plot(d, y, 'linewidth', 2);