function [delta_f, f0] = propagateErrorDxi(fun_handle, x0_vals, dx_vals)
% fun_handle should be a scalar function of x0. x0 may have multiple
% entries. dx_vals are the variations in each of the entries. delta_f is
% the expected variation in the function output 

% error propagation method
% df^2 ~ ((df/dx1)dx1)^2 + ((df/dx2)dx2)^2 + ...


% montecarlo methods

num_draws = 1000;

f_vals = zeros(num_draws, 1);
for ii = 1:num_draws
    x_prime = x0_vals + dx_vals.*randn(size(x0_vals));
    f_vals(ii) = fun_handle(x_prime);
end

delta_f = std(f_vals);
f0 = mean(f_vals);