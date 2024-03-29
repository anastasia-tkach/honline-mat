function [df] = my_gradient(f, x)

e = 10e-10;

%% If f is a scalar function

if isscalar(f(x))
    
    df = zeros(size(x'));
    
    for i = 1:length(x)
        delta_x = zeros(size(x));
        delta_x(i) = e;
        df(i) = (f(x + delta_x) - f(x - delta_x))/ (2 * e);
    end
    return;
end
    
%% If f is a vector function
    
if isvector(f(x))
    
    df = zeros(length(f(x)), length(x));
    
    for i = 1:length(x)
        delta_x = zeros(size(x));
        delta_x(i) = e;
        df(:, i) = (f(x + delta_x) - f(x - delta_x))/ (2 * e);
    end
    return;
end

%% If f is a matrix function
    
if ismatrix(f(x))
    
    df = zeros(size(f(x), 1), size(f(x), 2), length(x));
    
    for i = 1:length(x)
        delta_x = zeros(size(x));
        delta_x(i) = e;
        df(:, :, i) = (f(x + delta_x) - f(x - delta_x))/ (2 * e);
    end
    
end
