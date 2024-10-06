function [Y_gr, exp_X] = gr_cond(X_j, Y_j, theta_j)
%% [Y_gr, exp_X] = GR_COND(X_j, Y_j, theta_j)
%
% This function computes the regression of Y using conditional PDFs. This
% function is a sub-routine used for the estimation of the optimal
% smoothing parameter in the general regression approach.
%

%% Beginning

%% Compute size dataset
N_j = size(X_j, 1);

%% Pre-allocate space in memory
Y_gr     = zeros(N_j, 1);
exp_X    = zeros(N_j,1);
exp_X(1) = 1;

%% Assign first value
Y_gr(1) = Y_j(1);

%% Compute Y_hat using 'bisection' method
for jj = 2:N_j
  %% Compute mean distance between point being estimated and points in the conditional part
  tmppp     = gsubtract(X_j(1:jj-1,:),X_j(jj,:));
  dis_X     = mean(sum(tmppp.^2,2));
  exp_X(jj) = sum(exp(-2*sum(tmppp.^2,2)));
  
  %% Use GR to compute value of Y
  a = 0;
  b = 0;
  for m = 1:jj-1
    tmp = exp(-((X_j(jj,:) - X_j(m,:))*(X_j(jj,:) - X_j(m,:))')/(2*theta_j*dis_X+realmin));
    a   = a + Y_j(m)*tmp;
    b   = b + tmp;
  end
  
  Y_gr(jj) = a/(b+realmin);                 % Add 'realmin' in case of b = 0
end

end
%% END