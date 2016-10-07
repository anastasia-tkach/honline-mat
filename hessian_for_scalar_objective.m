function [h] = hessian_for_scalar_objective(F, J, H)

h = 2 * J' * J;
for i = 1:size(F, 1)
    h = h + 2 * F(i) * squeeze(H(i, :, :));
end