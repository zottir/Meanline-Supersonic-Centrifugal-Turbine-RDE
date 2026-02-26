function [unstarted] = checkUnstarting(x_shock,y_shock,Nb,r_in)
% Function to check the started condition, verifying that the bow shocks
% intersection point is placed meridionally after the LE radius.
% INPUT
% x_shock,y_shock: Bow shock coordinates [m];
% Nb: Number of blades;
% r_in: Inlet Radius
% OUTPUT
% unstarted: 0 if started condition, 1 if unstarted condition


% Copy the shock and rotate it
[x_shock_s1,y_shock_s1]= rotation(x_shock,y_shock,2*pi/Nb);

% Find Intersection
[x_i_shock,y_i_shock] = polyxpoly(x_shock_s1,y_shock_s1,x_shock,y_shock);

% Check the unstarting
if (x_i_shock^2 + y_i_shock^2) < r_in^2
    unstarted = 1;
else
    unstarted = 0;
end

end

