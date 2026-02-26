function [shock,blades,phi_shock,inters_points,refl_int,ds_bow,coll_shock] = bowShockLosses(M_in,gamma,R,inc_le,t_LE,ecc_LE,conc,alpha_in,alpha_out,Nb,r_in,r_out,chord,flagISO)
% Function to evaluate shock losses.
% INPUT
% M_in: Blade Inlet Mach Number [-];
% gamma: Specific Heat Ratio [-];
% R: Gas Constant [J/kgK];
% inc_le: Flow Incidence Angle [deg];
% t_LE: LE Thickness [m];
% ecc_LE: LE Eccentricity [-];
% conc: geom param. If equal to 1 deflect the flow in clock-wise direction,
%                   If equal to 2 in counter clock-wise direction.
% alpha_in,alpha_out: Inlet, Outlet Flow Angle (Relative if rotor) [rad];
% Nb: Number of Blades;
% r_in,r_out: Inlet, Outlet Radius [m];
% chord: Blade Chord [m];
% flagISO: if 1, ds = 0;
% OUTPUT
% shock: Matrix of Bow Shock Coordinates [x_shock; y_shock] [m];
% blades: Matrix with Approximate Blade Coordinates [x_blade; y_blade];
% phi_shock: Bow Shock Angle Distrinution [rad];
% inters_points: Matrix (2 x 2) containing upper and lower bow shock
% intersection points [x_i_upperBow,y_i_upperBow; x_i_lowerBow,y_i_lowerBow]
% refl_int
% 1 if bow shock refl on following/previous blade
% 2 if bow shock refl on same blade
% 0 if bow shock intersects circumference
% ds_bow: Entropy Production [J/kgK];
% coll_shock: flag equal to 1 if collective shock present;

%% Generate the approximate blade
[blades] = approximate_blade(alpha_in,alpha_out,r_in,r_out,chord,Nb,conc);
% Generate the leading edge
[x_LE_shape,y_LE_shape] = leadingEdge(0,t_LE/2,ecc_LE,gamma);

if flagISO == 1
    ds_bow = 0;
    shock = [];
    phi_shock = [];
    inters_points = [];
    refl_int = [];
    coll_shock = 0;
else

    %% Bow shock
    optionFinalSlopeCorrect = 1;
    % if 1, correct the slope of Moeckel shock based on corner analogy
    [deltax,x_shock_up,y_shock_up,x_shock_down,y_shock_down,phi_up,phi_down,~,SB,S] = bowed_shockWave(M_in,t_LE,x_LE_shape,y_LE_shape,alpha_in,inc_le,ecc_LE,gamma,R,6,1,optionFinalSlopeCorrect);

    % Traslation of the shock
    [x_shock_up,y_shock_up,~,~,~,~,~,~] = shockTraslation(S,SB,x_shock_up,y_shock_up,deltax,x_LE_shape,y_LE_shape,alpha_in,r_in);
    [x_shock_down,y_shock_down,~,~,~,~,~,~] = shockTraslation(S,SB,x_shock_down,y_shock_down,deltax,x_LE_shape,y_LE_shape,alpha_in,r_in);

    % Intersection
    [x_i_upperBow,y_i_upperBow,x_i_lowerBow,y_i_lowerBow,refl_up,refl_down] = intersection_bow_shock(x_shock_up,y_shock_up,x_shock_down,y_shock_down,blades,r_out);

    % Cut
    [x_shock_up,y_shock_up,phi_up] = shockCut(x_shock_up,y_shock_up,phi_up,x_i_upperBow,y_i_upperBow);
    [x_shock_down,y_shock_down,phi_down] = shockCut(x_shock_down,y_shock_down,phi_down,x_i_lowerBow,y_i_lowerBow);

    y_shock = [flip(y_shock_up) y_shock_down];
    x_shock = [flip(x_shock_up) x_shock_down];
    phi_shock = [flip(phi_up),phi_down];

    [ds_bow] = shock_Entropy(phi_shock,M_in,gamma,R,y_shock);

    % Outputs
    inters_points = [x_i_upperBow,y_i_upperBow; x_i_lowerBow,y_i_lowerBow];
    refl_int = [refl_up,refl_down];

    %% Check Unstarting
    [coll_shock] = checkUnstarting(x_shock,y_shock,Nb,r_in);

    shock = [x_shock;y_shock];

end
end