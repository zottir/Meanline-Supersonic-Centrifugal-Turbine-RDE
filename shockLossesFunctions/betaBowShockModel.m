function [beta_phiS] = betaBowShockModel(M1, phiS_ratio, gamma, R)
% Evaluation of the final bow shock slope with the simplyfied model based
% on the bow shock analogy with a corner oblique shock
%
% INPUT
% M1:               Inlet Mach Number
% gamma:            Specific Heat Ratio
% R:                Specific Gas Constant [J/kgK]
% phiS_ratio = 6;   Costant for the correlation to calculate theta oblique shock
%
% OUTPUT
% beta_phiS:        Oblique Shock Angle analogus to bow shock [rad]

%% Moeckel Theory
[phi_S,~,~,~,~,~] = Moeckel(M1, gamma, R, 1);

%% Correlation for theta, analogy corner oblique shock
theta_phiS = (90 - rad2deg(phi_S))/mean(phiS_ratio);

%% Oblique Shock Slope Calculation
beta_phiS = obliqueShockMachTheta(M1, deg2rad(theta_phiS), 1, gamma);

end

