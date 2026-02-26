function [blades] = approximate_blade(alpha_in,alpha_out,r_in,r_out,chord,Nb,conc)
% Function to approximately draw the blade starting from a cubic line,
% which respects the passage between 2 end points and the flow angles.
% INPUTS:
% alpha_in: Inlet Flow Angle [rad];
% alpha_out: Outlet Flow Angle [rad];
% r_in: Inlet Radius [m];
% r_out: Outlet Radius [m];
% chord: Blade Chord [m]
% Nb: Number of Blades;
% conc: default 1, 2 if rotor shape not good
% 
% OUTPUTS:
% blades: matrix (6 x n) with blades coordinates
%      | x      |
%      | y      |
%      | x_up   |
%      | y_up   |
%      | x_down |
%      | y_down |
% x: blade x coordinates;
% y; blade y coordinates;
% x_up: upper blade x coordinates;
% y_up: upper blade x coordinates;
% x_down: lower blade x coordinates;
% y_down: lower blade y coordinates;

%% Blade Design
phi = acos(0.5*(r_in^2+r_out^2-chord^2)/(r_in*r_out));
gamma = acos(0.5*(r_out^2-r_in^2+chord^2)/(chord*r_out));
sigma = asin(sin(phi)*r_out/chord);

if conc == 1

    % Write the polynomial in a local coordinate frame (x_loc,nu)
    A = chord*[1 chord; 2 3*chord];
    b = [-tan(sigma-alpha_in);tan(gamma-alpha_out)-tan(sigma-alpha_in)];
    pol = A\b;
    x_loc = linspace(0,chord,100);
    y_loc = tan(sigma-alpha_in)*x_loc + pol(1)*x_loc.^2 + pol(2)*x_loc.^3;

    % Conversion of the local to globa coordinate system
    R = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];
    sol = R(-sigma)*[x_loc;y_loc] + [r_in;0];

elseif conc == 2

    % Write the polynomial in a local coordinate frame (x_loc,nu)
    A = [chord chord^2;...
        2*chord 3*chord^2];
    B = [tan(sigma + alpha_in);tan(-gamma-alpha_out)+tan(sigma + alpha_in)];
    pol = A\B;
    x_loc = linspace(0,chord,100);
    y_loc = -tan(sigma + alpha_in)*x_loc + pol(1)*x_loc.^2 + pol(2)*x_loc.^3;

    % Conversion of the local to globa coordinate system
    R = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];
    sol = R(sigma)*[x_loc;y_loc] + [r_in;0];
    
end

x = sol(1,:);
y = sol(2,:);


%% Periodicity
dtheta = 2*pi/Nb;
sol_up = R(dtheta)*sol ;
x_up = sol_up(1,:);
y_up = sol_up(2,:);

sol_down = R(-dtheta)*sol ;
x_down = sol_down(1,:);
y_down = sol_down(2,:);

blades = [x; y; x_up; y_up; x_down; y_down];
end