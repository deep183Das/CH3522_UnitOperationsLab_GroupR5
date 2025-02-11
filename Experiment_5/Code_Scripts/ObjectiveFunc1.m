function [res] = ObjectiveFunc1(x,t,y)
N = length(y);
res = zeros(N,1);
for i = 1:N
    % residuals in y and y^hat : that's all which needs to be provided
    % as a vector.

    % for each sample compute it and pass it to the nonlinear function
    % It then itself will compute the sum squared error internally.
    res(i) = y(i) - (x(1) + (x(2)*erf((t(i)-x(3))/x(4))));
end