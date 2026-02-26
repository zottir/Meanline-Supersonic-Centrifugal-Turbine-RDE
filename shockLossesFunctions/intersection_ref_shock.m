function [x_ref_up,y_ref_up,x_ref_down,y_ref_down] = intersection_ref_shock(x_ref_up,y_ref_up,x_ref_down,y_ref_down,blades,x_i_upperBow,x_i_lowerBow,r_out)
% Function to evaluate the intersection point of the reflected shock and
% cut what is not needed
% INPUT
% x_ref_up,y_ref_up: upper reflected shock coordinates[m]
% x_ref_down,y_ref_down: : lower reflected shock coordinates [m];
% blades: matrix (6 x n) with blades coordinates
%      | x      |
%      | y      |
%      | x_up   |
%      | y_up   |
%      | x_down |
%      | y_down |
% x_i_upperBow,x_i_lowerBow: x coordinate of the up/low bow shock inters;
% r_out: Outlet radius [m]
%
% OUTPUT
% x_ref_up,y_ref_up: upper reflected shock coordinates [m];
% x_ref_down,y_ref_down: lower reflected shock coordinates [m];

%% Intersection Upper Reflected Shock
% Intersection with the lower blade
[x_i_upperRefl,y_i_upperRefl] = polyxpoly(x_ref_up,y_ref_up,blades(1,:),blades(2,:));
if ~isempty(x_i_upperRefl)
    x_i_upperRefl = x_i_upperRefl(1);
    y_i_upperRefl = y_i_upperRefl(1);
else
    % Intersection with the same blade
    [x_i_upperRefl,y_i_upperRefl] = polyxpoly(x_ref_up,y_ref_up,blades(3,:),blades(4,:));
    if ~isempty(x_i_upperRefl)
        % Don't consider the 1st inters if it's the bow shock incidence point
        condit = abs(x_i_upperRefl-x_i_upperBow) >  + 10^-5;
        x_i_upperRefl = x_i_upperRefl(condit);
        y_i_upperRefl = y_i_upperRefl(condit);
        % Take the first element of the remaining vector
        if ~isempty(x_i_upperRefl)
            x_i_upperRefl = x_i_upperRefl(1);
            y_i_upperRefl = y_i_upperRefl(1);
        else
            % Intersection with the circumference
            [x_i_upperRefl,y_i_upperRefl] = polyxpoly(x_ref_up,y_ref_up,r_out*cos(linspace(0,2*pi,500)),r_out*sin(linspace(0,2*pi,500)));
        end
    else
        % Intersection with the circumference
        [x_i_upperRefl,y_i_upperRefl] = polyxpoly(x_ref_up,y_ref_up,r_out*cos(linspace(0,2*pi,500)),r_out*sin(linspace(0,2*pi,500)));
    end
end

[~,ind] = min(sqrt((x_ref_up-x_i_upperRefl).^2 + (y_ref_up-y_i_upperRefl).^2));
x_ref_up = x_ref_up(1:ind);
y_ref_up = y_ref_up(1:ind);



%% Intersection Lower Reflected Shock
% Intersection with the upper blade
[x_i_lowerRefl,y_i_lowerRefl] = polyxpoly(x_ref_down,y_ref_down,blades(1,:),blades(2,:));
if ~isempty(x_i_lowerRefl)
    x_i_lowerRefl = x_i_lowerRefl(1);
    y_i_lowerRefl = y_i_lowerRefl(1);
else
    % Intersection with the same blade
    [x_i_lowerRefl,y_i_lowerRefl] = polyxpoly(x_ref_down,y_ref_down,blades(5,:),blades(5,:));
    if ~isempty(x_i_lowerRefl)
        % Don't consider the 1st inters if it's the bow shock incidence point
        condit = abs(x_i_lowerRefl-x_i_lowerBow) > 10^-5;
        x_i_lowerRefl = x_i_lowerRefl(condit);
        y_i_lowerRefl = y_i_lowerRefl(condit);
        % Take the first element of the remaining vector
        if ~isempty(x_i_upperRefl)
            x_i_lowerRefl = x_i_lowerRefl(1);
            y_i_lowerRefl = y_i_lowerRefl(1);
        else
            % Intersection with the circumference
            [x_i_lowerRefl,y_i_lowerRefl] = polyxpoly(x_ref_down,y_ref_down,r_out*cos(linspace(0,2*pi,500)),r_out*sin(linspace(0,2*pi,500)));
        end
    else
         % Intersection with the circumference
         [x_i_lowerRefl,y_i_lowerRefl] = polyxpoly(x_ref_down,y_ref_down,r_out*cos(linspace(0,2*pi,500)),r_out*sin(linspace(0,2*pi,500)));
    end
end

[~,ind] = min((x_ref_down-x_i_lowerRefl).^2 + (y_ref_down-y_i_lowerRefl).^2);
x_ref_down = x_ref_down(1:ind);
y_ref_down = y_ref_down(1:ind);
end



