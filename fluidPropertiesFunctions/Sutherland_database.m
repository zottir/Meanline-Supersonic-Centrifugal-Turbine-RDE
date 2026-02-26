function [mu, k] = Sutherland_database(T,gas_name, dati_GasViscosityDatabase)
% Function to calculate viscosity and themral conductivity based in
% Sutherland model
%
% INPUT:
% T: Reference temperature [K]
% gas_name: String with gas name
% dati_GasViscosityDatabase: table of the type: 
% | Fluid   mu0   T0mu   Smu   k0   T0k   Sk |
% |------------------------------------------|
% |------------------------------------------|
%
% OUTPUT:
% mu: viscosity [Pa*s]
% k: Thermal Conductivity [W/m*K]

%% Import data
B_GasViscosityDatabase=table2array(dati_GasViscosityDatabase(:,2:7));
gas_dataNames=convertCharsToStrings(table2array(dati_GasViscosityDatabase(1:end,1)))';
data_line=find(gas_dataNames==gas_name);
mu_0=B_GasViscosityDatabase(data_line, 1);
T_0mu=B_GasViscosityDatabase(data_line, 2);
S_mu=B_GasViscosityDatabase(data_line, 3);
k_0=B_GasViscosityDatabase(data_line, 4);
T_0k=B_GasViscosityDatabase(data_line, 5);
S_k=B_GasViscosityDatabase(data_line, 6);

%% Viscosity and Thermal Conductivity Evaluation
mu = mu_0 * ((T/T_0mu)^(3/2)) * ((T_0mu+S_mu)/(T+S_mu));
k = k_0 * ((T/T_0k)^(3/2)) * ((T_0k+S_k)/(T+S_k));
end

