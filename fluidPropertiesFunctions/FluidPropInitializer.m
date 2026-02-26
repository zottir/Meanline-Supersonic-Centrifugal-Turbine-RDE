 % clear;
 % close all;
 % clc;

%% Fluid Properties Initializer for Gamma, R, cp
% Input File for cp, gamma, R calculation
dati_flow = './fluidPropertiesFunctions/fluids/RDE_simplfied.xlsx';
dati=readtable(dati_flow);

% Import data
FlowData=table2array(dati(:,2:4));
molecule = convertCharsToStrings(table2array(dati(1:end,1)))';
MolecularMasses=FlowData(:,1);     % Molar Masses g/mol
MassFractions=FlowData(:,2);       % Mass Fractions
MolarFractions=FlowData(:,3);       % Molar Fractions

% Import of NASA constants from themodynamic data: given the molecule, find
% the 5 coefficients. Actually each coefficient is defined for 2 T ranges.
% const_up: T > 1000 K
% const_low: T < 1000 K

zconst_up=[];
zconst_low=[];
for zMol=1:max(size(molecule))
    molecola=convertStringsToChars(molecule(zMol));
    [zconst_up_mol,zconst_low_mol] = importMoleculeData(molecola);
    zconst_up=[zconst_up; zconst_up_mol];
    zconst_low=[zconst_low; zconst_low_mol];
end

%% Fluid Property Initializer for Viscosity
dati_GasViscosityDatabase = readtable('GasViscosityDatabase.xlsx');
ViscosityMixtureFormula = 2;  %1 Graham, 2 Wilke-Golubov, 3 Herning and Zipperer, 4 Brokaw, 5 Davidson
