% Author CLiu
% Date 07-03-19
% Computational simulation for validating the contrast model in LSCI 
% conventional model - Model1  --> mid-sized vessel but has been used
% everywhere 
% New derived model  - Model2 --> large vessel
%                    - Model3 --> parenchyma
clc
clear all
T=5000; % exposure time [us]   
tauc=1:1:50000;  %decorrelation time [us]
beta=1;  % optical parameter 
x=T./tauc;  
for i=1:length(tauc)
K1(i)=Model1(beta,tauc(i),T);
K2(i)=Model2(beta,tauc(i),T);
K3(i)=Model3(beta,tauc(i),T);
Ksimp1(i)=Simp1(beta,tauc(i),T);
end
plot(tauc,K1,'g');
hold on
plot(tauc,K2,'r-');
plot(tauc,K3,'b-');
%% Error caused by wrong g1 model------------------------------
T = 5;
x1 = 0;
x2 = 10^10;
beta = 1;
%large vessels
taucBase =0.05;  
trueRF = 0.1:0.1:3;
taucResp = taucBase ./ trueRF;
Kbase = Model2(beta,taucBase,T);
for i=1:length(taucResp)
    Kresp(i) = Model2(beta,taucResp(i),T);
end
taucBase_m = taucModel1(Kbase,beta,x1,x2,T);
for i=1:length(Kresp)
    taucResp_m(i) = taucModel1(Kresp(i),beta,x1,x2,T);
end
RF_large_n1 = taucBase_m ./ taucResp_m;
RElarge_n1 = ( RF_large_m-trueRF)./trueRF ;  

RF_large_simp= (Kbase.^2)./ (Kresp.^2);
RE_large_simp= (RF_large_simp - trueRF)./trueRF;

%parenchyma region
taucBase =1;  
trueRF = 0.1:0.1:3;
taucResp = taucBase ./ trueRF;
Kbase = Model3(beta,taucBase,T);
for i=1:length(taucResp)
    Kresp(i) = Model3(beta,taucResp(i),T);
end

taucBase_m = taucModel1(Kbase,beta,x1,x2,T);
for i=1:length(Kresp)
    taucResp_m(i) = taucModel1(Kresp(i),beta,x1,x2,T);
end

RF_par_n1 = taucBase_m ./ taucResp_m;
RE_par_n1 = ( RF_par_n1-trueRF)./trueRF ;  
RF_par_simp= (Kbase.^2)./ (Kresp.^2);
RE_par_simp= (RF_par_simp - trueRF)./trueRF;

%mid-sized vessels
taucBase =0.2;  
trueRF = 0.1:0.1:3;
taucResp = taucBase ./ trueRF;
Kbase = Model1(beta,taucBase,T);
for i=1:length(taucResp)
    Kresp(i) = Model1(beta,taucResp(i),T);
end
RF_mid_simp= (Kbase.^2)./ (Kresp.^2);
RE_mid_simp= (RF_mid_simp - trueRF)./trueRF;

%% Error caused by beta
% method: trueRF --correct model (Model1 Model2 Model3,correct beta) -- contrast K---
% correct g1 model wrong beta(beta=1) ---measured RF
% the baseline tauc is in ms units
% large vessel : taucBase = 0.05 ms
% small vessel : taucBase = 0.2 ms
% Parenchyma   : taucBase = 1 ms
% remember to use the correct Model for each situation
clc
clear all
T = 5;
x1 = 0;
x2 = 10^10;
taucBase = 1;
%taucBase = 0.2;
%taucBase = 0.05;
trueRF = 0.1:0.01:3;
taucResp = taucBase ./ trueRF;
truebeta =0.1:0.1:1;

for i=1:length(truebeta)
KBase = Model1(truebeta(i), taucBase, T);
taucBase_m(i) = taucModel1(KBase,1,x1,x2,T);
end
tic
for i=1:length(truebeta)
    i
for j = 1:length(taucResp)
    KResp = Model1(truebeta(i), taucResp(j),T);
    taucResp_m(i,j) = taucModel1(KResp,1,x1,x2,T);
end
end
toc
for i=1:length(truebeta)
RF_m(i,:) = taucBase_m(i) ./ taucResp_m(i,:);
RE(i,:) = (RF_m(i,:) - trueRF )./ trueRF;
end

RF_large_beta1=RF_m;
RFm_par_beta1=RF_m;
RFm_mid_beta1=RF_m;

RE_large_beta = RE;
RE_mid_beta= RE;
RE_par_beta = RE;
%% Error caused by static scattering 
% method: trueRF --correct model (Model1 Model2 Model3,correct beta,correct rho) -- 
%contrast K---correct g1 model (beta=1), rho=1 ---measured RF
% the baseline tauc is in ms units
% large vessel : taucBase = 0.05 ms
% small vessel : taucBase = 0.2 ms
% Parenchyma   : taucBase = 1 ms
% make sure to use correct model for each case
clc
clear all
T=5;
beta=1;
x1=0;
x2=10^10;
taucBase =0.05; taucBase = 1; taucBase = 0.2  ;
trueRF = 0.1:0.01:3;
taucResp = taucBase ./ trueRF;
truerho=0.1:0.1:1;

for i=1:length(truerho)
    i
    KBaseCom = SModel3( beta,truerho(i),T,taucBase );
    taucB_m(i) = taucModel3( KBaseCom,beta,x1,x2,T );
    for j=1:length(taucResp)
        kRespCom = SModel3(beta,truerho(i),T,taucResp(j));   
        taucResp_m(i,j)=taucModel3( kRespCom,beta,x1,x2,T );
    end
end
for i=1:length(truerho)
RF_m(i,:)= taucB_m(i) ./ taucResp_m(i,:);
RE(i,:)= (RF_m(i,:)-trueRF) ./ trueRF; 
end
RE_par_rho1=RE;
%RE_mid_rho1=RE;
%RE_large_rho1=RE;













