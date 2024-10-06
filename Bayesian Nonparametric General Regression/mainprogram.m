%% Beginning

%% Close open windows, clear memory and screen
close all; clear; clc;

%% Set seed for random number generator (for reproducibility)
%     seed_rng = rng(0, 'twister');
%%

% Input data (design matrix)
                         
y1  = xlsread('data','sheet1','a2:a102');
y2  = xlsread('data','sheet1','b2:b102');
y3  = xlsread('data','sheet1','c2:c102');

N   = size(y1,1);                 % Size dataset
X   = [y1 y2 y3]; 
d   = size(X,2);                  % Number of design variables

% Output data (Real data)
y_n = xlsread('data','sheet1','d2:d102'); 

%% Normalize data (mean=0, std=1)

% Normalize design matrix
mean_x = mean(X);
std_x  = std(X);
X_norm = zeros(N,d);
for i = 1:d
  X_norm(:,i) = (X(:,i) - mean_x(i))./std_x(i);
end

% Normalize measurements
mean_y   = mean(y_n);
std_y    = std(y_n);
Y_n_norm = (y_n - mean_y)./std_y;

%% Define training and testing datasets
N_tr = floor(0.7*N);          % Size training dataset (70% of the data)取比它小的整数
N_te = N - N_tr;              % Size testing dataset (30% of the data)

idx_te         = (1:N)';
idx_tr         = (sort(randperm(N, N_tr)))';    % Indices training set
idx_te(idx_tr) = [];                            % Indices testing set

% Training set
X_tr        = X(idx_tr,:);
X_norm_tr   = X_norm(idx_tr,:);
y_tr        = y_n(idx_tr);
y_n_tr      = y_n(idx_tr);
Y_n_norm_tr = Y_n_norm(idx_tr);

% Testing set
X_te        = X(idx_te,:);
X_norm_te   = X_norm(idx_te,:);
y_te        = y_n(idx_te);
y_n_te      = y_n(idx_te);
Y_n_norm_te = Y_n_norm(idx_te);

%% Model matrix
model      = dec2bin(0:2^(d)-1)-'0';
model(1,:) = [];
N_mod      = size(model, 1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%   Training stage   %%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters prior PDF
l_b = 0;                            % Lower boundary
u_b = 50;                           % Upper boundary

%% Pre-allocate space in memory
theta_hat   = zeros(2,N_mod);       % Optimal parameters for each model class
Y_n_norm_gr = zeros(N_tr,N_mod);    % Regression of Y computed via GR for each model class (Normalized values)
Y_hat       = zeros(N_tr,N_mod);    % Regression of Y computed via GR for each model class
MSE         = zeros(N_mod, 1);      % Mean-square error
like_k      = zeros(N_mod, 1);      % Likelihood function value for each model class
P_Y_X_C     = zeros(N_mod, 1);      % Evidence of each model class

%% Options minimization routine
options = optimset('Display','off','TolX',1e-4,'TolFun',1e-4);

%% Print information on screen
fprintf('Model    Theta_1     Theta_2    P(Y|theta,X,C)   MSE        p(Y|X,C)\n');
fprintf('-------------------------------------------------------------------------\n');
tic;

%% Main code
for k = 1:N_mod
  % Define handle function used to compute optimal smoothing parameter
  % (General regression)
  neg_log_pos = @(theta) min_neg_log_pos(l_b, u_b, ...
                                            X_norm_tr(:,model(k,:)==1), ...
                                                            Y_n_norm_tr, theta);

  % Compute optimal value of variance of the kernel used in GR using
  % minimization routine
  while theta_hat(1,k) <= 0
      theta_hat(1,k)=fminsearch(neg_log_pos, 0.1*rand, options);
  end

  % Compute optimal value of prediction-error variance
  [Y_tmp, exp_X] = gr_cond(X_norm_tr(:,model(k,:)==1), Y_n_norm_tr, theta_hat(1,k));
  MSE_tmp        = (Y_n_norm_tr - Y_tmp)'*diag(exp_X)*(Y_n_norm_tr - Y_tmp);
  theta_hat(2,k) = MSE_tmp / N_tr;

  % Compute likelihood function value
  like_k(k) = exp(log_like_fun(Y_n_norm_tr, Y_tmp, theta_hat(2,k), exp_X));
%   like_k(k) = (1/((2*pi*theta_hat(2,k))^(0.5*N_tr))) * ...
%                                               exp(-MSE_tmp/(2*theta_hat(2,k)));

  % Compute evidence of model class M_k
  P_Y_X_C(k) = evid_class(l_b, u_b, X_norm_tr(:,model(k,:)==1), ...
                                                  Y_n_norm_tr, theta_hat(:,k));

  % Compute regression of Y in X using the whole dataset
  Y_n_norm_gr(:,k) = gr_equation(X_norm_tr(:,model(k,:)==1), Y_n_norm_tr, ...
                                    X_norm_tr(:,model(k,:)==1), Y_n_norm_tr, ...
                                                            theta_hat(1,k), 0);

  % Re-scale regression of Y in X
  Y_hat(:,k) = (Y_n_norm_gr(:,k).*std_y) + mean_y;

  % Compute Mean-squared error
  MSE(k) = (y_tr - Y_hat(:,k))'*(y_tr - Y_hat(:,k))/N_tr;

  % Update information on screen
  fprintf('%2d       %2.4f      %2.4f     %2.4e       %2.4f     %2.4e\n', ...
              k, theta_hat(1,k), theta_hat(2,k), like_k(k), MSE(k), P_Y_X_C(k));
end
fprintf('-------------------------------------------------------------------------\n');
fprintf('\n');
total_time = toc;

%% Compute plausilibty of each model class
P_C_Y_X = (P_Y_X_C ./ N_mod) ./ sum(P_Y_X_C ./ N_mod);

%% Identify best model
best_model = find(P_C_Y_X == max(P_C_Y_X));
sigma2_gr  = theta_hat(1,best_model);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%   Testing stage   %%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Identify design variables best model
X_te_best = X_norm_te(:,model(best_model,:)==1);
X_tr_best = X_norm_tr(:,model(best_model,:)==1);

%% Compute regression of Y_te using optimal smoothing parameter computed in training stage
Y_n_norm_te_gr = gr_equation(X_te_best, Y_n_norm_te, ...
                                    X_tr_best, Y_n_norm_tr, sigma2_gr, 1);

%% Re-scale output data
Y_te_gr = (Y_n_norm_te_gr .* std_y) + mean_y;
