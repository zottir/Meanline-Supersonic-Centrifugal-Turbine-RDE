function [phi_S, phi_C, Ps_by_P0_c, C, sigma_shock, lambda_S] = Moeckel(M_1, gamma, R, geom)
% Function to calculate Moeckel Parameters to draw the bow shock
% INPUT:
% M_1: INlet Mach Number [];
% gamma: Specifi Heat Capaity Ratio [];
% R: Specifi Gas Constant [J/kgK]
% geom: 1 if 2D geometry, 2 if axial-symmetric
%
% OUTPUT
% phi_S: Shock Angle at sonic point [rad] 
% phi_C: Shock Angle in Sonic Line Centroid [rad] 
% Ps_by_p0_c: Total Pressure Jump in Centroid [-] 
% C: Moeckel Constant
% sigma_shock: Isentropic Section Variation to have Sonic Conditions
% lambda_S: Flow Angle after Sonic Point [rad]
%
% Notation (Moeckel vs Anderson)
% phi - beta : Shock Angles
% lambda - theta:  Flow Angle

if geom == 1  
     % 2D Geometry 
    yc_by_yS = 1/2;
elseif geom == 2    
    % Axial Symmetric Geometry
    yc_by_yS = 2/3;
end

Cp = R*gamma/(gamma-1);
    
%% Iterative Cycle

M_2_obj = 1;
beta = asin(1/M_1);
err = 1;
toll = 1e-6;
Nmax = 1000;
k = 1;

while err>toll && k<Nmax  
    M_1n = M_1*sin(beta);
    M_2n = sqrt((M_1n^2 + (2/(gamma-1))) / ((((2*gamma)/(gamma - 1))*(M_1n^2)) - 1));
    theta = atan((2/tan(beta)) * ((M_1n^2 - 1) / ((M_1^2)*(gamma + cos(2*beta)) + 2)));
    M_2 = (M_2n)/(sin(beta-theta));
    
    err = abs((M_2 - M_2_obj)/(M_2_obj));
    beta = beta + 0.05*((M_2 - M_2_obj)/(M_2_obj));
    k = k+1;   
end

phi_S = beta;
lambda_S = theta;

%% Centroid Angle
beta_shock = sqrt(M_1^2 - 1);
phi_C = abs(atan(sqrt((((beta_shock*tan(phi_S))^2) + (yc_by_yS^2) - 1)/((beta_shock*yc_by_yS)^2))));
C = beta_shock*((beta_shock*tan(phi_S)) - sqrt(((beta_shock*tan(phi_S))^2) - 1));
sigma_shock = 1/(sqrt((1./(M_1).^2).*((((2)./(gamma+1)) .* (1 + ((gamma - 1)./(2)).*M_1.^2)).^((gamma+1)./(gamma-1)))));

%% Centroid Total Pressure Jump Ratio
%Normal shock relations
M_1n_C = M_1 * sin(phi_C);
p2_by_p1 = 1 + (((2*gamma)/(gamma+1))*(M_1n_C.^2 - 1));
rho2_by_rho1 = ((gamma+1).*M_1n_C.^2)./((2 + ((gamma-1)*M_1n_C.^2)));
T2_by_T1 = p2_by_p1/rho2_by_rho1;

%Entropy variation accross shock e total pressure drop
deltas = Cp*log(T2_by_T1) - R*log(p2_by_p1);
pt2_by_pt1 = exp((-deltas)./(R));

Ps_by_P0_c = pt2_by_pt1;

end

