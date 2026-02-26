function [XX, YY] = leadingEdge(beta_in_deg,radius_1,ecc_LE, gamma)
% Function to draw the leading edge of the blade
% INPUT
% beta_in_deg: Inlet Flow Angle [deg]
% radius_1: LE thickness radius [m];
% ecc_LE: LE eccentricity;
% gamma: Specific Heat Capacity Ratio []
%
% OUTPUT
% XX, YY: x and y coordinates of the LE

%% Fluid Data
a=((gamma+1)./(gamma-1));

%% Draw Parameters Evaluation
beta_in_rad=deg2rad(beta_in_deg);
radius= radius_1/cos(beta_in_rad); %ottengo il distanziamento assiale dei due profili

%% Elliptic edges

M=[0; radius_1];
N=[0; -radius_1];

C_LE=[(M(1,1)+N(1,1))/2; (M(2,1)+N(2,1))/2];
b_LE=(M(2,1)-N(2,1))/2;
a_LE=sqrt((b_LE^2)/(1-ecc_LE^2));

% Leading Edge

x_LE=linspace(C_LE(1,1)-a_LE, C_LE(1,1),500);
y_LE_up=C_LE(2,1)+sqrt((1 - (x_LE-C_LE(1,1)).^2./(a_LE.^2))*(b_LE^2));
y_LE_down=C_LE(2,1)-sqrt((1 - (x_LE-C_LE(1,1)).^2./(a_LE.^2))*(b_LE^2));
x_LE_upDown=[x_LE(end:-1:1), x_LE];
y_LE_upDown=real([y_LE_up(end:-1:1), y_LE_down]);

beta_in_rot=-beta_in_rad;
mat_rot=[cos(-beta_in_rot)   -sin(-beta_in_rot);
        sin(-beta_in_rot)    cos(-beta_in_rot)];
    
xy_LE_rot=mat_rot*[x_LE_upDown; y_LE_upDown];
x_LE_upDown=xy_LE_rot(1,:);
y_LE_upDown=xy_LE_rot(2,:);

[~,I] = max(y_LE_upDown);
y_LE_upDown = y_LE_upDown - (y_LE_upDown(I));
x_LE_upDown = x_LE_upDown - (x_LE_upDown(I));

% %% Horizontal Section
% 
% x0_SS=x_LE_upDown(1);
% y0_SS=y_LE_upDown(1);
% 
% x0_PS=x_LE_upDown(end);
% y0_PS=y_LE_upDown(end);
% 
% x_SS=linspace(x0_SS, x0_SS+25*radius);
% y_SS=tan(beta_in_rad)*(x_SS-x0_SS) + y0_SS;
% 
% x_PS=linspace(x0_PS, x0_PS+25*radius);
% y_PS=tan(beta_in_rad)*(x_PS-x0_PS) + y0_PS;
% 
% %% Points for the profile (with horiz section)
% 
% XX= [x_SS x_LE_upDown flip(x_PS)]';        %rounding_out_stat(1,end:-1:1)
% YY= [y_SS y_LE_upDown flip(y_PS)]';
XX= x_LE_upDown';        %rounding_out_stat(1,end:-1:1)
YY= y_LE_upDown';
[min_XX, I_min_XX]=min(XX);
XX= XX - min_XX;
YY= YY - YY(I_min_XX);

end

