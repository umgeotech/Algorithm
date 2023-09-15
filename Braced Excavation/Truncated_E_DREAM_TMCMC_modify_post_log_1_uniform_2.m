% Sampling algorithm
function [x,ln_S]=Truncated_E_DREAM_TMCMC_modify_post_log_1_uniform_2 (log_like_fun,N,log_mean,log_deviation,nnn,low_bound,up_bound,exponential);

%% initialing

p = 0; ln_S = 0;
%% drawing prior samples
for i=1:2;
x(i,:) =low_bound(i)+(up_bound(i)-low_bound(i))*rand(1,N); 
end
i=3;
x(i,:)=lognrnd(log_mean,log_deviation,[1,N]);
for i = 1:N,
    log_like(i,1) = feval(log_like_fun,x(:,i),exponential);
end
[delta,c,c_star,Pm,p_g]=deal(3,0.1,1e-12,0.9,0.2);
lambda=unifrnd(-c,c,N,1);
%% compute the weights
% adaptively choose p
while p<1,
    low_p = p; up_p = 2; old_p = p;
    while up_p-low_p>1e-6
        current_p = (low_p + up_p)/2;
        temp = exp((current_p-p)*(log_like-max(log_like)));
        cov_temp = std(temp)/mean(temp);
        if cov_temp > 1
            up_p = current_p;
        else
            low_p = current_p;
        end
    end
    p = current_p;
    if p > 1, break; end % breakout if approaching the final stage
    weight = temp/sum(temp); % weights are normalized
    ln_S = ln_S+log(mean(temp))+(p-old_p)*max(log_like);

    [~,posi_max]=max(weight);
%% do MCMC
 sam_ind = deterministicR((1:N),weight); current_x = x; current_log_like = log_like;
 for i=1:N,
now_ind = sam_ind(i);
[kk] = randperm(N);
ii = kk(1:2);a = ii(1); b = ii(2);
z=rand(1,nnn);
A=find(z<Pm); 
d_star=numel(A); 
if d_star==0,[~,A]=min(z);d_star=1;end 
gamma_d=2.38/sqrt(2*delta*d_star);
g=randsample([gamma_d 1],1,'true',[1-p_g p_g]); 
dX(1:nnn,i)=0;
dX(A,i)=(1+lambda(i))*(g*(current_x(A,posi_max)-current_x(A,now_ind))+g*(current_x(A,a)-current_x(A,b)))+c_star*randn(d_star,1);
x_c=current_x(:,now_ind)+dX(:,i);
for j=1:2;
if x_c(j)>up_bound(j)
x_c(j)=low_bound(j)+(up_bound(j)-low_bound(j))*rand(1);
end
if x_c(j)<low_bound(j)
x_c(j)=low_bound(j)+(up_bound(j)-low_bound(j))*rand(1);
end
end

old_x=current_x(:,now_ind);
j=3;
current_prior=lognpdf(old_x(j),log_mean,log_deviation);
x_prior=lognpdf(x_c(j),log_mean,log_deviation);
        log_like_c = feval(log_like_fun,x_c,exponential);
        r = exp(p*(log_like_c - current_log_like(now_ind))+log(x_prior)-log(current_prior));

        if r > rand,
            x(:,i) = x_c; current_x(:,now_ind) = x_c; current_log_like(now_ind) = log_like_c; log_like(i) = log_like_c;
        else
            x(:,i) = current_x(:,now_ind); log_like(i) = current_log_like(now_ind);
        end
    end
end
%% finalize the final stage
temp = exp((1-old_p)*(log_like-max(log_like)));
weight = temp/sum(temp);
[~,posi_max]=max(weight);
ln_S = ln_S+log(mean(temp))+(1-old_p)*max(log_like);
sam_ind = deterministicR((1:N),weight); current_x = x; current_log_like = log_like;
for i=1:N,
now_ind = sam_ind(i);
[kk] = randperm(N);
ii = kk(1:2);a = ii(1); b = ii(2);
z=rand(1,nnn);
A=find(z<Pm); 
d_star=numel(A); 
if d_star==0,[~,A]=min(z);d_star=1;end 
gamma_d=2.38/sqrt(2*delta*d_star);
g=randsample([gamma_d 1],1,'true',[1-p_g p_g]); 
dX(1:nnn,i)=0;
dX(A,i)=(1+lambda(i))*(g*(current_x(A,posi_max)-current_x(A,now_ind))+g*(current_x(A,a)-current_x(A,b)))+c_star*randn(d_star,1);
x_c=current_x(:,now_ind)+dX(:,i);
for j=1:2;
if x_c(j)>up_bound(j)
x_c(j)=low_bound(j)+(up_bound(j)-low_bound(j))*rand(1);
end
if x_c(j)<low_bound(j)
x_c(j)=low_bound(j)+(up_bound(j)-low_bound(j))*rand(1);
end
end

old_x=current_x(:,now_ind);
j=3;
current_prior=lognpdf(old_x(j),log_mean,log_deviation);
x_prior=lognpdf(x_c(j),log_mean,log_deviation);
        log_like_c = feval(log_like_fun,x_c,exponential);
        r = exp(log_like_c - current_log_like(now_ind)+log(x_prior)-log(current_prior));
    if r > rand, 
        x(:,i) = x_c; current_x(:,now_ind) = x_c; current_log_like(now_ind) = log_like_c; log_like(i) = log_like_c;
    else
        x(:,i) = current_x(:,now_ind); log_like(i) = current_log_like(now_ind);
    end
end
