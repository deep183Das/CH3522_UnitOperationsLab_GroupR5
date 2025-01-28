function [res] = ConcFunc(x,y,t)
N = length(t);
res = zeros(N,1);
for i = 1:N
    % residuals in y and y^hat : that's all which needs to be provided
    % as a vector

    % for each sample compute it and pass it to the nonlinear function
    % It then itself will compute the sum squared error internally
    res(i) = y(i) - (x(1) + x(2)*exp(-x(3)*t(i)));
end