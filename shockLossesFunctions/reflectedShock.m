function [delta,x_ref,y_ref,phi2,M1,detached,no_refl_Mn] = reflectedShock(M_in,alpha_in,phi_in,blades,x_i_Bow,y_i_Bow,refl,r_out,gamma,pos)
% Function to draw the reflected shock
% INPUTS
% M0: Inlet Mach Number [-];
% alpha0: Inlet Flow Angle [rad];
% phi0: Bow Shock Angle Distribution [rad];
% blades: matrix (6 x n) with blades coordinates
%      | x      |
%      | y      |
%      | x_up   |
%      | y_up   |
%      | x_down |
%      | y_down |
% x_i_Bow,y_i_Bow: Intersection point coordinates of the bow shock;
% refl: flag equal to 0 if no reflection, 1 if bow shock intersects the
% following/previous blade, 2 if bow shock intersects the same blade;
% r_out: Outlet Radius [-];
% gamma: Specific Heat Ratio [-];
% pos: 1 for upper reflection, 2 for lower reflection
% All distributions from upper branch to lower one
% OUTPUT
% delta: Oblique Shock Angle with respect to the radial direction [rad]
% x_ref,y_ref: coordinates of the oblique shock points [m]
% M1: Mach number after the bow shock
% detached: flag equal to 1 if the reflected shock is detached
% no_refl_Mn: flag equal to 1 if the reflected shock is not present thanks
% to Mout(bow_shock) < 1.

%% Reflected Shock upper branch
if pos == 1
    detached = 0;
    no_refl_Mn = 0;
    if refl ~=0
        % 0: before bowed shock
        % 1: after bowed shock - before reflected shock
        % 2: after reflected shock

        % Post-shock properties
        Mn_in = M_in*sin(phi_in(1));
        theta1 = atan(2*(Mn_in^2 - 1)/(tan(phi_in(1)) * (M_in^2*(gamma + cos(2*phi_in(1))) + 2)));
        Mn1 = sqrt((2+(gamma-1)*Mn_in^2)/(2*gamma*Mn_in^2 - (gamma-1)));
        M1 = Mn1/sin(phi_in(1)-theta1);

        % Find the slope of the intersected blade
        if refl == 1
            [~,ind] = min(abs(blades(3,:) - x_i_Bow));
            beta = atan((blades(4,ind+1) - blades(4,ind))/(blades(3,ind+1) - blades(3,ind)));
        else
            [~,ind] = min(abs(blades(1,:) - x_i_Bow));
            beta = atan((blades(2,ind+1) - blades(2,ind))/(blades(1,ind+1) - blades(1,ind)));
        end
        theta2 = (theta1 - alpha_in) - beta;

        options = optimoptions('fsolve', 'Display', 'off');
        fphi = @(phi) 2*((M1*sin(phi))^2 - 1)/(tan(phi) * (M1^2*(gamma + cos(2*phi)) + 2)) - tan(abs(theta2));
        [phi2, ~, exitflag,~] = fsolve(fphi,phi_in(1),options);

        % Check a solution exists
        if (exitflag <= 0) || (M1*sin(phi2)<1)
            if exitflag <= 0
            detached = 1;
            end
            if M1*sin(phi2)<1
            no_refl_Mn = 1;
            end
            phi2 = [];
            delta = [];
            x_ref = [];
            y_ref = [];
            M1 = [];
        else
            delta = theta2 - phi2 + beta;
            x_ref = linspace(x_i_Bow,r_out*1.05,500);
            y_ref = y_i_Bow + delta*(x_ref - x_i_Bow);
        end
    else
        phi2 = [];
        delta = [];
        x_ref = [];
        y_ref = [];
        M1 = [];
    end
end



%% Reflected Shock lower branch
if pos == 2
    detached = 0;
    no_refl_Mn = 0;
    if refl ~=0
        % 0: before bowed shock
        % 1: after bowed shock - before reflected shock
        % 2: after reflected shock

        % Post-shock properties
        Mn_in = M_in*sin(phi_in(end));

        theta1 = atan(2*(Mn_in^2 - 1)/(tan(phi_in(end)) * (M_in^2*(gamma + cos(2*phi_in(end))) + 2)));
        Mn1 = sqrt((2+(gamma-1)*Mn_in^2)/(2*gamma*Mn_in^2 - (gamma-1)));
        M1 = Mn1/sin(phi_in(end)-theta1);
        
        % Find the slope of the intersected blade
        if refl == 1
            [~,ind] = min(abs(blades(5,:) - x_i_Bow));
            beta = atan((blades(6,ind+1) - blades(6,ind))/(blades(5,ind+1) - blades(5,ind)));
        else
            [~,ind] = min(abs(blades(1,:) - x_i_Bow));
            beta = atan((blades(2,ind+1) - blades(2,ind))/(blades(1,ind+1) - blades(1,ind)));
        end
        theta2 = (theta1 + alpha_in) + beta;

        options = optimoptions('fsolve', 'Display', 'off');
        fphi = @(phi) 2*((M1*sin(phi))^2 - 1)/(tan(phi) * (M1^2*(gamma + cos(2*phi)) + 2)) - tan(abs(theta2));
        [phi2, ~, exitflag,~] = fsolve(fphi,phi_in(end),options);

        % Check a solution exists
        if (exitflag <= 0) || (M1*sin(phi2)<1)
            if exitflag <= 0
                detached = 1;
            end
            if M1*sin(phi2)<1
                no_refl_Mn = 1;
            end
            delta = [];
            x_ref = [];
            y_ref = [];
            phi2 = [];
            M1 = [];
        else
            delta = phi2 - abs(theta1 + alpha_in);
            x_ref = linspace(x_i_Bow,r_out*1.05,500);
            y_ref = y_i_Bow + delta*(x_ref - x_i_Bow);
        end
    else
        delta = [];
        x_ref = [];
        y_ref = [];
        phi2 = [];
        M1 = [];
    end
end
end