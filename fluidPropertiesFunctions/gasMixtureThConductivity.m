function [k_gasMixture] = gasMixtureThConductivity(T,dati_mixture, dati_GasViscosityDatabase)
% Function to calculate mixture thermal conductivity with Mason-Saxena
% model
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
%
% OUTPUT
% k_gasMixture: Mixture thermal conductivity


%% Import data gas mixture
B=table2array(dati_mixture(:,2:5));
molecole=convertCharsToStrings(table2array(dati_mixture(1:end,1)))';
Mm=B(:,1)';      
X=B(:,4)';       

%% Calculus molecule thermal conductivity with Sutherland
k_vec=[];
mu_vec=[];
clear i
for i=1:max(size(molecole))
    gasName_i=molecole(i);
    [mu_i, k_i] = Sutherland_database(T,gasName_i, dati_GasViscosityDatabase);
    k_vec=[k_vec k_i];
    mu_vec=[mu_vec mu_i];   
end

%% Mason-Saxena Model(1958)
clear i j
G=[];
for i=1:max(size(molecole)) 
    for j=1:max(size(molecole))    
        G_part=(1 + ((((mu_vec(i)*Mm(j))/(mu_vec(j)*Mm(i)))^(1/2))*((Mm(i)/Mm(j))^(1/4))))^(2);
        G(i,j)=(1.065/(2*sqrt(2)))*((1 + (Mm(i)/Mm(j)))^(-1/2))*G_part;    
    end  
end
G=G-diag(diag(G));  % Make null diagonal element to have i != j
k_gasMixture=sum((k_vec)./(1 + (sum((((X'*ones(1,max(size(molecole))))'./(X'*ones(1,max(size(molecole))))).*G), 2)')));

end

