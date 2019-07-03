% Contrast model for parenchyma flow measurement 
% get the tau_c from contrast using complete model g1(tau)=exp(-(tau/tau_c)^{0.5});

% Inputs:
%      K      -Contrast value obtained from data collection
%      Beta   -beta value
%      x0    - initial guess in fminsearch function
%      T      -exposure time (ms)
%      rho    - dynamic component 

%  Output:
%      tauc   -decorrelation time 

%  Ratio=x=T/tau_c

%------------- BEGIN CODE --------------

function [ tauc ] =taucSModel3( K,beta,x0,T,rho )

syms x 

   fun = @(x) abs( sqrt( beta*( rho^2* (4*x.*exp(-2*sqrt(x))+6*sqrt(x).*exp(-2*sqrt(x)) +2*x+3*exp(-2*sqrt(x))-3 )./(2*x.^2) +...
       8*rho*(1-rho)* (2*x.*exp(-sqrt(x))+6*sqrt(x).*exp(-sqrt(x))+x+6*exp(-sqrt(x))- 6  )./x.^2  +(1-rho)^2 ))  -K);
  
  x=fminsearch(fun,x0);
  tauc=T./x;
  


end

