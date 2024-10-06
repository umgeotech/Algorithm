function log_like = log_like_fun(y, y_hat, v2, dis_X)
%% log_like = log_like_fun(y, y_hat, v2, dis_X)
%
% This function is used to compute the log-likelihood function.

%% Beginning

%% Compute size dataset
N = length(y);

%% Initialize value of log_like. Take into account that p(y_1 | y_0, v1, v2, X, C) = 1
log_like = 0;

%% Compute log_like for each conditional PDF and add the values
for i = 1:N
  %% Acummulate value of logarithm of likelihood function
  log_like = log_like - 0.5*log(2*pi*(v2/dis_X(i))) - ...
                                      (((y(i) - y_hat(i))^2)/(2*(v2/dis_X(i))));
  
end

end
%% END