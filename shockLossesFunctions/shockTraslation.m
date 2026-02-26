function [x_shock,y_shock,XX_prof,YY_prof,x_sb,y_sb,x_s,y_s] = shockTraslation(S,SB,x_shock,y_shock,deltax,XX_prof,YY_prof,alpha_in,r_in)
% Function to transformate the local shock reference frame into the global
% one.
% INPUT
% SB: Body Sonic Point Coordinates
% |x_sb_down x_sb_up |
% |y_sb_down y_sb_up |
% S: Sonic Point Coordinates
% |x_s_down x_s_up |
% |y_s_down y_s_up |
% x_shock, y_shock: Shock coordinates in the shock reference frame
% deltax: Distance between the shock foremost point and blade LE
% XX_prof,YY_prof: LE elliptical/circumf profiles
% alpha_in: Inlet Flow Angle [rad];
% r_in: Inlet radius
%
%OUTPUT
% Same quantities in the global reference frame


%% Translation of the shock in (0,0)
SB(1,:) = SB(1,:) - min(x_shock) - deltax;
S(1,:) = S(1,:) - min(x_shock) - deltax;
x_shock = x_shock  - min(x_shock) - deltax;
XX_prof = XX_prof - min(XX_prof);

% Rotation of the shock to match the angle of attack
[XX_prof,YY_prof] = rotation(XX_prof', YY_prof', -alpha_in);
[x_shock,y_shock] = rotation(x_shock, y_shock, -alpha_in);
[x_sb,y_sb] = rotation(SB(1,:),SB(2,:), -alpha_in);
[x_s,y_s] = rotation(S(1,:),S(2,:), -alpha_in);

% Translation of the points in the good radius
x_shock = x_shock + r_in;
x_sb = x_sb + r_in;
x_s = x_s + r_in;
XX_prof = XX_prof + r_in;

end