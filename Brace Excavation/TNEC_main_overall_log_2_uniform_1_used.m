% Main code for calculations
% Take the case of updating three unknown parameters using the measured deflection data in the third stage as an example.
clc;
close all;
clear;
tic;
times_number=10;
N=5000;
burnIn=1;
COV=0.3;
nnn=3;
exponential=4;
for times=1:times_number;
log_like_fun = 'TNEC_post_3'; 

low_bound=[0,0];
up_bound=[20,2e4];
normal_mean=11.7; 
normal_deviation=COV*normal_mean;

M=normal_mean;
V=normal_deviation.^2;
log_mean=log(M.^2./sqrt(V+M.^2));
log_deviation=sqrt(log(V./M.^2+1));

[x,ln_S]=Truncated_E_DREAM_TMCMC_modify_post_log_1_uniform_2 (log_like_fun,N,log_mean,log_deviation,nnn,low_bound,up_bound,exponential);
mu1=mean(x(1,burnIn:N));
S1= std (x(1,burnIn:N));
mu2=mean(x(2,burnIn:N));
S2= std (x(2,burnIn:N));
mu3=mean(x(3,burnIn:N));
S3= std (x(3,burnIn:N));
TNEC_results_post(1,2*times-1:2*times)=[mu1,S1];
TNEC_results_post(2,2*times-1:2*times)=[mu2,S2];
TNEC_results_post(3,2*times-1:2*times)=[mu3,S3];
LNS_post(1,2*times-1:2*times)=ln_S;

x_error(times,:)=x(1,:);
x_cf(times,:)=x(2,:);
x_a(times,:)=x(3,:);
end
toc;   

% plot
min=0;
max=20;
space=0.02;
uuu=min:space:max;
nBins=length(uuu);
sampleBins=linspace(min,max,nBins);
figure;
counts1= hist (x(1,burnIn:N), sampleBins);
bar(sampleBins, counts1/space/sum(counts1), 'k');
xlabel('TMCMC samples' ); ylabel( 'Posterior function');

min=0;
max=10000;
space=5;
uuu=min:space:max;
nBins=length(uuu);
sampleBins=linspace(min,max,nBins);
figure;
counts2= hist (x(2,burnIn:N), sampleBins);
bar(sampleBins, counts2/space/sum(counts2), 'k');
xlabel('TMCMC samples' ); ylabel( 'Posterior function');

min=0;
max=20;
space=0.1;
uuu=min:space:max;
nBins=length(uuu);
sampleBins=linspace(min,max,nBins);
figure;
counts3= hist (x(3,burnIn:N), sampleBins);
bar(sampleBins, counts3/space/sum(counts3), 'k');
xlabel('TMCMC samples' ); ylabel( 'Posterior function');