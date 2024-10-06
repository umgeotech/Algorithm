function P_D_M = evid_class(l_b, u_b, X, Y, Theta)
%% P_D_M = EVID_CLASS(l_b, u_b, X, Y, Theta)
%
% This function computes the evidence of a model class defined by X using
% asymptotic expansion.
%
%% Beginning

%% Handle function for negative log-posterior PDF
target = @(theta) -log(sigma2_mar_pdf(l_b, u_b, X, Y, theta));

%% Compute Hessian matrix
H = hessian(target, Theta(1));

%% Compute evidence of model class M_k
P_D_M = exp(-1*target(Theta(1)))*sqrt(2*pi/det(H));

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%   Definition of functions used in this script   %%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Marginal PDF of variance of the kernel used in the General Regression approach
function p_sigma = sigma2_mar_pdf(l, u, X_j, Y_j, t)
%% p_sigma = SIGMA2_MAR_PDF(l, u, X_j, Y_j, t)
%
% This function computes the marginal PDF of the variance of the kernel used in
% the General Regression approach.
%
% (*) INPUT DATA:
%     - l:   Lower boundary prior PDF (real number)
%     - u:   Upper boundary prior PDF (real number)
%     - X_j: Design variables (Normalized values) (N_j x d matrix)
%     - Y_j: Measurements (Normalized values) (N_j x 1 vector)
%     - t:   Variance of the Gaussian kernel (positive real number)
%
% (*) OUTPUT DATA:
%     - p_sigma: marginal PDF of the variance of the kernel used in the General
%                Regression approach (real number)
%


%% Beginning
N_j = size(X_j,1);            % Size dataset

%% Compute Y_hat
[Y_hat, exp_X] = gr_cond(X_j, Y_j, t);

%% Compute marginal PDF
C       = (2/((u-l)^2))*(1/(pi^(N_j/2)))*gamma(0.5*N_j + 1);
l_p_s   = log(C) + sum(0.5*log(exp_X)) - ...
                    (0.5*N_j + 1)*log(((Y_j-Y_hat)'*diag(exp_X)*(Y_j-Y_hat)));
p_sigma = exp(l_p_s);
end

%% Hessian matrix using finite difference method
function Hess = hessian(fun,x)
%% Hess = HESSIAN(fun,x)
%
% This function computes the numeric Hessian matrix of a multivariate
% function. The Hessian matrix is computed by finite difference method.
%
% (*) INPUT DATA:
%
% - fun: Handle function
% - x:   Point at which the hessian will be evaluated (dim x 1 vector)
%
% (*) OUTPUT DATA:
%
% - Hess: Numeric Hessian matrix (dim x dim matrix)
%


%% Beginning

%% Size of Hessian matrix
x   = x(:);
dim = length(x);

%% Create perturbation vectors:
delta_Theta = diag(0.01.*abs(x));

%% Evaluate the function at x
fx = feval(fun,x);diag(0.01.*abs(x))

%% Compute Hessian by central difference
Hess = zeros(dim);
for m = 1:dim
  %% Diagonal elements
  Hess(m,m) = (feval(fun, x + delta_Theta(:,m)) - ...
                              2*fx + feval(fun, x - delta_Theta(:,m))) ./ ...
                                          (delta_Theta(:,m)'*delta_Theta(:,m));

  if Hess(m,m) <= 0
    Hess(m,m) = eps;
  end

  %% Off-diagonal elements, use central difference method
  for n = m+1:dim
    delta_1   = delta_Theta(:,m);
    delta_2   = delta_Theta(:,n);
    Hess(m,n) = (feval(fun, x + delta_1 + delta_2) - ...
                      feval(fun, x + delta_1 - delta_2) - ...
                            feval(fun, x - delta_1 + delta_2) + ...
                                  feval(fun, x - delta_1 - delta_2)) ./ ...
                                      (4*(delta_Theta(m,m)*delta_Theta(n,n)));

    % Hess is a symmetric matrix
    Hess(n,m) = Hess(m,n);
  end
end

end

%% END