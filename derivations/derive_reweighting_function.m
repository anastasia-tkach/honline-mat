
d = 0:0.00001:20;
weight1 = zeros(length(d), 1);
weight0 = zeros(length(d), 1);

factor = 3.5;

for i = 1:length(d)
    w = (d(i) + 1e-3)^(-0.5);
    if (d(i) > 1e-3)       
        weight1(i) = w * factor;
    end
    weight0(i) = w * factor;
end
weight2 = weight1;
weight2(weight2 > 8) = 8;

figure; hold on;
plot(d, d, 'lineWidth', 2);
plot(d, weight0' .* d.^2, 'lineWidth', 2);
plot(d, weight1' .* d.^2, 'lineWidth', 2);
plot(d, weight2' .* d.^2, 'lineWidth', 2);
%legend('w = 3.5/sqrt(d + 1e-3)', 'if d < 1e-3, w = 0', 'if w > 10, w = 10');