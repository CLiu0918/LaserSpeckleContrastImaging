%Conventional contrast model 
%get the tau_c from contrast using complete model g1(tau)=exp(-tau/tau_c)
% Inputs:
%      K      -Contrast value obtained from data collection
%      Beta   -beta value
%      x0    - initial guess in fminsearch function
%      T      -exposure time (ms)
%      rho    - dynamic component 

%  Output:
%      tauc   -decorrelation time 

%------------- BEGIN CODE --------------
function [ tauc ] =taucSModel1( K,beta,x0,T,rho )

syms x 

  fun=@(x) abs(  sqrt( beta .* ( rho^2*( exp(-2*x) - 1 + 2*x)./(2*x.^2) +4*rho*(1-rho).*...
    (exp(-x)-1+x)./x.^2+(1-rho)^2) )-K  );
  
  x=fminsearch(fun,x0);
  tauc=T./x;  % ms
  
end

