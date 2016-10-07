function [segments] = shape_2D_autodiff(segments, beta)

segments{2}.local(2, 4) = beta(1);
segments{3}.local(2, 4) = beta(2);
segments{4}.local(2, 4) = beta(3);


