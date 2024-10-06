function Y_hat = gr_equation(X, Y, X_k, Y_k, theta, type)
%% Y_hat = GR_EQUATION(X, Y, X_k, Y_k, theta, type)
%
% This function computes the regression of Y in X using the sample values
% (X_k, Y_k).
%
%     - theta: Variance of the kernel used in General Regression
%              (positive real scalar)
%      - type: Type of problem:
%              * type = 0 ---> X == X_k
%              * type = 1 ---> X ~= X_k
%

%% Beginning

%% Compute size datasets
N   = size(Y,1);
N_k = size(X_k,1);

%% Compute distances
dis = zeros(N, 1);

for i = 1:N
  for j = 1:N_k
    dis(i) = dis(i) + (X(i,:) - X_k(j,:)) * (X(i,:) - X_k(j,:))';
  end
end

%% Pre-allocate space in memory
Y_hat = zeros(N,1);

%% Use GR equation to compute Y_hat
for j = 1:N
  a = 0;
  b = 0;
  sigma_1 = theta*dis(j)/N_k;
  for k = 1:N_k
    if type == 0
      if k ~= j
        D_2 = (X(j,:) - X_k(k,:)) * (X(j,:) - X_k(k,:))';
        a = a + Y_k(k)*exp(-D_2/(2*sigma_1));
        b = b + exp(-D_2/(2*sigma_1));
      end
    else
      D_2 = (X(j,:) - X_k(k,:)) * (X(j,:) - X_k(k,:))';
      a = a + Y_k(k)*exp(-D_2./(2*sigma_1));
      b = b + exp(-D_2./(2*sigma_1));
    end
  end
  Y_hat(j) = a/(b+realmin);
end
end
%% END