
d = 0:0.0001:10;
weight = zeros(length(d), 1);

for i = 1:length(d)
if (d(i) > 1e-3)
    w = 1/sqrt(d(i) + 1e-3);
    weight(i) = w * 3.5;
end
end

figure; plot(d, weight, 'lineWidth', 2); 