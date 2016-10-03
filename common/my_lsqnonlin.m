function [X, J] = my_lsqnonlin(function_handle, X, num_iters)

E_history = zeros(num_iters, 1);
lambda = 1;
E_previous = Inf;
for iter = 1:num_iters + 1

    [F, J] = function_handle(X);
    E = F' * F;

    if E <= E_previous
        
        %lambda = lambda / 2;
        lambda = lambda / 1.5;
        E_previous = E;
        F_previous = F;
        J_previous = J;
        X_previous = X;    
       
        E_history(iter) = E;
    else
        %lambda = lambda * 10;
        lambda = lambda * 5;
        E = E_previous;
        F = F_previous;
        J = J_previous;
        X = X_previous;

        E_history(iter) = E_history(iter - 1);
    end

    if iter == num_iters, break; end
    
    delta = - (J' * J + lambda * eye(size(J, 2), size(J, 2))) \ (J' * F);  
    X = X + delta;
end

%figure; plot(1:num_iters, E_history, 'lineWidth', 2);

