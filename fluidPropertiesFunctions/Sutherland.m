function [mu, k] = Sutherland(T,dati_mixture,dati_GasViscosityDatabase, ViscosityMixtureFormula)
% Function to calculate the constants of Sutherland’s formula through a 
% least-squares interpolation, given the viscosity and thermal conductivity
% of the mixture as temperature varies.
%
% INPUT 
% T: Temperature at which evaluate the properties [K]
% dati_mixture: table of the type:
% | Molecule   Molecular_Mass   Mass_Fraction   Molar_Fraction|
% |-----------------------------------------------------------|
% |-----------------------------------------------------------|
% dati_GasViscosityDatabase: table of the type: 
% | Fluid   mu0   T0mu   Smu   k0   T0k   Sk |
% |------------------------------------------|
% |------------------------------------------|
% ViscosityMixtureFormula: 1 Graham, 2 Wilke-Golubov, 3 Herning and Zipperer, 4 Brokaw, 5 Davidson
% 
% OUTPUT
% mu: Viscosity [Pa*s]
% k: Thermal Conductivity [W/m*K]

%% Evaluation of viscosity(T)
T_vec=linspace(350, 3500, 100);
mu_mixture_vec=[];
k_mixture_vec=[];
for i=1:max(size(T_vec))
    T=T_vec(i);
    [mu_gasMixture] = gasMixtureViscosity(T,dati_mixture, dati_GasViscosityDatabase, ViscosityMixtureFormula);
    [k_gasMixture] = gasMixtureThConductivity(T,dati_mixture, dati_GasViscosityDatabase);
    mu_mixture_vec=[mu_mixture_vec mu_gasMixture];
    k_mixture_vec=[k_mixture_vec k_gasMixture];
end

%% Evaluation lsq constants for viscosity curve fitting 
xdata=T_vec;
ydata=mu_mixture_vec*10^(5);    %moltiplico per 10^5 perché i valori sono piccolo e le tolleranza da imporre diventano ancora più piccole
F=@(x,xdata)x(1).*((xdata./x(2)).^(3/2)).*((x(2)+x(3))./(xdata+x(3))); 
x0=[1.716; 273; 111];

options = optimoptions('lsqcurvefit','FunctionTolerance', 10^(-20),'Display','none');
lb = [];
ub = [];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,xdata,ydata,lb,ub, options);

mu_0=x(1)*10^(-5);          %moltiplico per 10^5 per tornare ai valori assoluti corretti
T_0mu=x(2);
S_mu=x(3);

%% Evaluation lsq constants for thermal conductivity curve fitting  
xdata=T_vec;
ydata=k_mixture_vec;
F=@(x,xdata)x(1).*((xdata./x(2)).^(3/2)).*((x(2)+x(3))./(xdata+x(3))); 
x0=[0.0241; 273; 194];

options = optimoptions('lsqcurvefit','FunctionTolerance', 10^(-20),'Display','none');
lb = [];
ub = [];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,xdata,ydata,lb,ub, options);

k_0=x(1);
T_0k=x(2);
S_k=x(3);

mu = mu_0 * ((T/T_0mu)^(3/2)) * ((T_0mu+S_mu)/(T+S_mu));
k = k_0 * ((T/T_0k)^(3/2)) * ((T_0k+S_k)/(T+S_k));


end

