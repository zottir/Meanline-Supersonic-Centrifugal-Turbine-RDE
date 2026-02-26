function [ds_shock] = shock_Entropy(phi,M_in,gamma,R,y_shock)
% Function to calculate entropy across a bowed or oblique shock.
% INPUT:
% phi: Shock angle distribution [rad];
% M_in: Shock Inlet Mach Number [-];
% gamma: Specific Heat Ratio;
% R: Gas Constant [J/kgK];
% y_shock: y Shock Coordiates [m]
% OUTPUT
% ds_shock: Entropy production across the shock [J/kgK]

% Entropy Calculation
if ~isempty(phi)
cp = gamma*R/(gamma-1);
Mn_in = M_in*sin(phi);
dp = 1 + 2*gamma/(gamma+1)*(Mn_in.^2-1);
drho = (gamma+1)*Mn_in.^2./(2+(gamma-1)*Mn_in.^2);
dT = dp./drho;
ds_distribution = cp*log(dT) - R*log(dp);

% Integration along the channel
ds_shock = trapz(y_shock,ds_distribution)/(y_shock(end)-y_shock(1));
else
    ds_shock = 0;
end

end