function [refl_shock_up,refl_shock_down,ds_refl_shock,detached,no_refl_Mn] = reflectShockLosses(M_in,gamma,R,alpha_in,r_out,phi_shock,ds_bowed,blades,inters_points,refl_int,flagISO)
% Function to evaluate shock losses.
% INPUT
% M_in: Blade Inlet Mach Number [-];
% gamma: Specific Heat Ratio [-];
% R: Gas Constant [J/kgK];
% alpha_in: Inlet Flow Angle (Relative if rotor) [rad];
% r_out: Outlet Radius [m];
% phi_shock: Bow Shock Angle Distrinution [rad];
% blades: matrix (6 x n) with blades coordinates
%      | x      |
%      | y      |
%      | x_up   |
%      | y_up   |
%      | x_down |
%      | y_down |
% inters_points: Matrix (2 x 2) containing upper and lower bow shock
% intersection points [x_i_upperBow,y_i_upperBow; x_i_lowerBow,y_i_lowerBow]
% refl_int:
% 1 if bow shock refl on following/previous blade
% 2 if bow shock refl on same blade
% 0 if bow shock intersects circumference
% ds_shock: Entropy Production [J/kgK];
% flagISO: if 1, ds = 0;
%
% OUTPUT
% refl_shock_up: Matrix of Upper Reflected Shock Coordinates [x_refl_shock; y_refl_shock] [m];
% refl_shock_down: Matrix of Lower Reflected Shock Coordinates [x_refl_shock; y_refl_shock] [m];
% ds_refl_shock: Entropy Production [ds_refl_up, ds_refl_down] [J/kgK];
% detached: flag equal to 1 if the Reflected Shock is detached [detached_up,detached_down];
% no_refl_Mn: flag equal to 1 if the reflected shock is not present thanks
% to Mout(bow_shock) < 1.

if flagISO == 1
    ds_refl_down = 0;
    ds_refl_up = 0;
    x_refl_up = [];
    y_refl_up = [];
    x_refl_down = [];
    y_refl_down = [];
    detached_up = 0;
    detached_down = 0;
    no_refl_Mn_up = 0;
    no_refl_Mn_down= 0;
else
    %% Reflected shock
    refl_up = refl_int(1);
    refl_down = refl_int(2);
    x_i_upperBow = inters_points(1,1);
    y_i_upperBow = inters_points(1,2);
    x_i_lowerBow = inters_points(2,1);
    y_i_lowerBow = inters_points(2,2);

    % Draw the reflected shock
    [~,x_refl_up,y_refl_up,phi_refl_up,M1_up,detached_up,no_refl_Mn_up] = reflectedShock(M_in,alpha_in,phi_shock,blades,x_i_upperBow,y_i_upperBow,refl_up,r_out,gamma,1);
    [~,x_refl_down,y_refl_down,phi_refl_down,M1_down,detached_down,no_refl_Mn_down] = reflectedShock(M_in,alpha_in,phi_shock,blades,x_i_lowerBow,y_i_lowerBow,refl_down,r_out,gamma,2);

    % Intersection and cut
    [x_refl_up,y_refl_up,x_refl_down,y_refl_down] =...
        intersection_ref_shock(x_refl_up,y_refl_up,x_refl_down,y_refl_down,blades,x_i_upperBow,x_i_lowerBow,r_out);

    % Entropy Production
    if detached_up == 0
        phi_refl = phi_refl_up.*ones(size(y_refl_up));
        [ds_refl_up] = shock_Entropy(phi_refl,M1_up,gamma,R,y_refl_up);
    else
        % Assume the reflected detached shock produces the same amount of loss
        % of the relative bow shock branch (conservative assumption).
        ds_refl_up = 0.5*ds_bowed;
    end

    if detached_down == 0
        phi_refl = phi_refl_down.*ones(size(y_refl_down));
        [ds_refl_down] = shock_Entropy(phi_refl,M1_down,gamma,R,y_refl_down);
    else
        % Assume the reflected detached shock produces the same amount of loss
        % of the relative bow shock branch (conservative assumption).
        ds_refl_down = 0.5*ds_bowed;
    end
end

    ds_refl_shock = [ds_refl_up,ds_refl_down];
    refl_shock_up = [x_refl_up;y_refl_up];
    refl_shock_down = [x_refl_down;y_refl_down];
    detached = [detached_up,detached_down];
    no_refl_Mn = [no_refl_Mn_up,no_refl_Mn_down];

end