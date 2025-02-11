function [res] = ObjectiveFunc(x,X,Y,nu,ratio_of_Co_to_Ci)
N = length(Y);
res = zeros(N,1);
for i = 1:N
    % residuals in y and y^hat : that's all which needs to be provided
    % as a vector.

    % for each sample compute it and pass it to the nonlinear function
    % It then itself will compute the sum squared error internally.
    res(i) = ratio_of_Co_to_Ci(i) - 0.5*(1 + erf((((3*Y(i)*x(2))/(2*X))-1)/(2*sqrt((nu*x(2))/(X*x(1))))));
    % The model expression in the objective function that is being used in
    % Rosen model to analyze adsorption breakthrough in fixed bed system.
end