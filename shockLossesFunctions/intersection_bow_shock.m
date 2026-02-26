function [x_i_upperBow,y_i_upperBow,x_i_lowerBow,y_i_lowerBow,refl_up,refl_down] = intersection_bow_shock(x_shock_up, y_shock_up,x_shock_down, y_shock_down,blades,r_out)
% Function to find the first intersection between the bow shock wave and
% the blade/circumference
% INPUT
% x_shock_up, y_shock_up: Upper bow shock coordinates
% x_shock_down, y_shock_down: Lower bow shock coordinates
% blades: matrix (6 x n) with blades coordinates
%      | x      |  
%      | y      |
%      | x_up   |
%      | y_up   | 
%      | x_down |
%      | y_down |
% r_out: Outlet radius of the blade
%
% OUTPUT
% x_i_upperBow,y_i_upperBow: Upper bow shock intersection coordinates
% x_i_lowerBow,y_i_lowerBow: Lower bow shock intersection coordinates
% unstarted: 1 if unstarted, 0 if not
% refl: 
% 1 if bow shock refl on following/previous blade
% 2 if bow shock refl on same blade
% 0 if bow shock intersects circumference
        
%% Intersection Upper Bow Shock
% 1) Intersect the upper blade;
[x_i_upperBow, y_i_upperBow] = polyxpoly(x_shock_up, y_shock_up, blades(3,:), blades(4,:));
if ~isempty(x_i_upperBow)
    % Take only the first Intersection
    x_i_upperBow = x_i_upperBow(1);
    y_i_upperBow = y_i_upperBow(1);
    refl_up = 1;
else
    % 2) Intersect the same blade;
    [x_i_upperBow, y_i_upperBow] = polyxpoly(x_shock_up, y_shock_up, blades(1,:), blades(2,:));
    if ~isempty(x_i_upperBow)
        % Take only the first Intersection
        x_i_upperBow = x_i_upperBow(1);
        y_i_upperBow = y_i_upperBow(1);
        refl_up = 2;
    else
        % 3) Intersect the circumference
        [x_i_upperBow, y_i_upperBow] = polyxpoly(x_shock_up, y_shock_up,r_out*cos(linspace(0,2*pi,600)),r_out*sin(linspace(0,2*pi,600)));
        if ~isempty(x_i_upperBow)
        % Take only the first Intersection
        x_i_upperBow = x_i_upperBow(1);
        y_i_upperBow = y_i_upperBow(1);
        else
            x_i_upperBow = x_shock_up(end);
            y_i_upperBow = y_shock_up(end);
        end
        refl_up = 0;
    end
end



%% Intersection Lower Bow Shock
% 1) Intersect the lower blade;
[x_i_lowerBow, y_i_lowerBow] = polyxpoly(x_shock_down, y_shock_down, blades(5,:), blades(6,:));
if ~isempty(x_i_lowerBow)
    % Take only the first Intersection
    x_i_lowerBow = x_i_lowerBow(1);
    y_i_lowerBow = y_i_lowerBow(1);
    refl_down = 1;
else
    % 2) Intersect the same blade;
    [x_i_lowerBow, y_i_lowerBow] = polyxpoly(x_shock_down, y_shock_down, blades(1,:), blades(2,:));
    if ~isempty(x_i_lowerBow)
        % Take only the first Intersection
        x_i_lowerBow = x_i_lowerBow(1);
        y_i_lowerBow = y_i_lowerBow(1);
        refl_down = 2;
    else
        % 3) Intersect the circumference
        [x_i_lowerBow, y_i_lowerBow] = polyxpoly(x_shock_down, y_shock_down,r_out*cos(linspace(0,2*pi,500)),r_out*sin(linspace(0,2*pi,500)));
        if ~isempty(x_i_lowerBow)
        % Take only the first Intersection
        x_i_lowerBow = x_i_lowerBow(1);
        y_i_lowerBow = y_i_lowerBow(1);    
        else
            x_i_lowerBow = x_shock_down(end);
            y_i_lowerBow = y_shock_down(end);    
        end
        refl_down = 0;
    end
end

end