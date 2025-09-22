%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%%     Triple Oxygen Isotope Parameter estimation by MCMC  %
%%                          EQUISETUM
%                   By Vincent Hare, October 2025                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uses MCMC toolbox (mjlaine.github.io/mcmcstat/), and
% prediction_plot_from_n.m, equisetum.mat, ssfun.m,
% model_equisetum3_plot.m, prediction_plot_from_n.m, model_equisetum3.m,
% and - optionally - my_priorfun.m, a custom prior function. 
%
% Estimated parameters are as follows:
% [qk ea es T d18Omw m]; 
%
% Dataset "NAME.mat" has the stucture:
%   NAME.ydata : d18O D'17Obw     
%   NAME.xdata : [1:1:length(NAME.mat)]     
%   NAME.sigma  : ud18O uD'17Obw (1 sigma)    

clear all

load equisetum.mat
par = [1.0350; 0.511; 1.0; 1.5; 1; 25.0; -7.0; 20];
nsimu = 1000000;
method      = 'dram';
adaptint    = 50;
N = length(data.xdata);
p = 6;
mse = ssfun(par,data)/(N-p);

params = {
%      name,    init,min, max, mu, sig, target?, local?
{'ak', par(1), 0.5, 1.5, par(1), 0.0001,   0,      1}
{'qk', par(2), 0.4, 6, par(2), 0.005,   1,      1}
{'ea', par(3), 0.001, 3.0, par(3), 0.5, 1, 1}
{'es', par(4), 0.001, 3.0, par(4), 0.5, 1, 1}
{'y', par(5), 0.1, 5.0, par(5), 0.1, 0, 1}
{'T', par(6), 0.001, 50, par(6), 1.0, 1, 1}
{'d18Omw', par(7), -15, 0, par(7), 1.0, 1, 1}
%{'m', par(8), -1, 1, par(8), 0.01, 1, 1}
{'n', par(8), N, 150, par(8), Inf, 1, 1}
};

model.ssfun=@ssfun;
model.sigma2 = 1; 

options.updatesigma = 1;
options.method      = method;
options.nsimu       = nsimu; 

[res,chain,s2chain] = mcmcrun(model,data,params,options);
chainstats(chain,res)

%Calculate AIC/BIC
fixed_ak = par(1);
fixed_y  = par(5);
par_mean = mean(chain); 
par_mean = [
    fixed_ak;
    par_mean(1);  % qk
    par_mean(2);  % ea
    par_mean(3);  % es
    fixed_y;
    par_mean(4);  % T
    par_mean(5);  % d18Omw
    par_mean(6)   % m
];
ss = ssfun(par_mean,data);
s2 = mean(s2chain(length(s2chain)-50:length(s2chain)));
logL = -0.5 * N * log(2 * pi) - 0.5 * N * log(s2) - 0.5 * ss / s2;
AIC = -2 * logL + 2 * p;
BIC = -2 * logL + p * log(N);

fprintf('AIC: %.2f\n', AIC);
fprintf('BIC: %.2f\n', BIC);

figure(1)
mcmcplot(chain,[],res,'chainpanel');
figure(2)
mcmcplot(chain,[],res,'denspanel',2);
figure(3)
mcmcplot(chain,[],res,'pairs');
prediction_plot_from_n(round(mean(chain(:,6))) ,chain, data, @model_equisetum3_plot, 10000);