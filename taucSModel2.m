% Contrast model for large vessel flow measurement
% get the tau_c from contrast using complete model g1(tau)=exp(-(tau/tau_c)^2);

% Inputs:
%      K      -Contrast value obtained from data collection
%      Beta   -beta value
%      x1,x2  -boundary conditions for function "fminbnd"
%      T      -exposure time (ms)
%      rho    - dynamic component 

%  Output:
%      tauc   -decorrelation time 

%  Ratio=x=T/tau_c

%------------- BEGIN CODE -------------- 
function [ tauc ] =taucSModel2(K,beta,x1,x2,T,rho )

syms x 

 
  fun=@(x) abs ( sqrt( beta* (rho^2 * ( (exp(-2*x.^2)+sqrt(2*pi).*x.*erf(sqrt(2)*x)-1)./(2*x.^2) )+...
      2*rho*(1-rho)* (exp(-x.^2)+sqrt(pi)*x.*erf(x)-1)  /x.^2  +(1-rho)^2  )  ) -K);
  
  Ratio=fminbnd(fun,x1,x2);
  tauc=T./Ratio;
  


end

