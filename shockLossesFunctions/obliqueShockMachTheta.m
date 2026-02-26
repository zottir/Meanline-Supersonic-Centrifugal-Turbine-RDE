function [beta] = obliqueShockMachTheta(M,theta, delta, gamma)
% Function to evaluate oblique shock properties
% INPUT
% M: Inlet Mach Number []
% theta: Corner Angle [rad] 
% delta:  1 if weak solution,0 if strong solution
% gamma: Specific Heat CApacity Ratio []
%
% OUTPUT
%beta: Oblique Shock Angle [rad] 

lamda_sup = (((M^2 - 1)^2) - (3*(1 + ((gamma-1)/2)*(M^2))*(1 + ((gamma+1)/2)*(M^2))*((tan(theta))^2)))^(1/2);
X_sup = (((M^2 - 1)^3) - 9*(1 + ((gamma-1)/2)*(M^2))*(1 + ((gamma-1)/2)*(M^2) + ((gamma+1)/4)*(M^4))*((tan(theta))^2))/(lamda_sup^3);
beta = atan(((M^2) - 1 + (2*lamda_sup*cos((4*pi*delta + ((acos(X_sup))))/3)))/(3*(1 + ((gamma-1)/2)*(M^2))*tan(theta)));   

end
