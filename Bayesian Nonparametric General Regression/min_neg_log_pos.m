function neg_log_pos = min_neg_log_pos(l_b, u_b, X, Y, theta)
%% neg_log_pos = MIN_NEG_LOG_POS(l_b, u_b, X, Y, theta)

%% Compute size dataset
N = size(X,1);

if (theta < l_b) || (theta > u_b)
  %% Compute negative log-posterior PDF
  neg_log_pos = realmax;
else
  %% Compute regression of Y using conditional PDFs and General Regression
  [Y_hat, exp_X] = gr_cond(X, Y, theta);

  %% Compute optimal value of the prediction-error variance
  R2       = (Y - Y_hat)'*diag(exp_X)*(Y - Y_hat);
  sigma2_2 = R2 / N;
  
  %% Compute negative log-posterior PDF
  neg_log_pos = 2*log(u_b-l_b) - log_like_fun(Y, Y_hat, sigma2_2, exp_X);
%   neg_log_pos = 2*log(u_b-l_b) + 0.5*N*log(2*pi*sigma2_2) + (R2/(2*sigma2_2));
end

%% In case of divergence, assign large positive value
if isnan(neg_log_pos) || isinf(neg_log_pos)
  neg_log_pos = realmax;
end

end

%% END