
clear;
close all;
clc;


load('32_57.mat')
writeHubShroud = 1;

%% STATOR CFD Setup
dr_s = 0.5*(r1-r0);
ri_s = r0 - dr_s;
ro_s = r1 + 1.5*c_stator;

T_fun = @(T) T0t - 0.5/cp*((r0*rho0*V0/ri_s)*R*(T/p0).*(T0./T).^(gamma/(gamma-1))).^2 - T;

% Serach of the Guess value for fsolve
Ti_vec_s = linspace(1.1*T0,2*T0,100);
T_fun_s = T_fun(Ti_vec_s);
sign_change_idx = find(T_fun_s(1:end-1).*T_fun_s(2:end) < 0);
Ti_zero_s = zeros(size(sign_change_idx));
for ki = 1:length(sign_change_idx)
    i = sign_change_idx(ki);
    % Valori ai due lati del cambio di segno
    x1 = Ti_vec_s(i);
    x2 = Ti_vec_s(i+1);
    y1 = T_fun_s(i);
    y2 = T_fun_s(i+1);
    % Linear interpolation
    Ti_zero_s(ki) = x1 - y1*(x2 - x1)/(y2 - y1);
end
%Choice of the supersonic solution
Ti_s = Ti_zero_s(1);
clear("Ti_zero_s","sign_change_idx","ki","x1",'"x2',"y1","y2")

% Stator BC
Ti_s = fzero(T_fun,Ti_s);
Vi_s = sqrt(2*cp*(T0t-Ti_s));
Mi_s = Vi_s/sqrt(gamma*R*Ti_s);
pi_s = p0*(Ti_s/T0)^(gamma/(gamma-1));
rhoi_s = pi_s/(R*Ti_s);



%% ROTOR CFD Setup
dr_r = 0.2*(r3-r2);
ri_r = r2 - dr_r;
ro_r = r3 + 1.5*c_rotor;

T_funR = @(T) T2t - 0.5/cp*( (V2t*r2/ri_r)^2 + ( (m*R*T/(2*pi*ri_r*b2*p2)) .* (T2./T).^(gamma/(gamma-1))  ).^2 ) - T;

% Serach of the Guess value for fsolve
Ti_vec_r = linspace(1*T2,2*T2,100);
T_fun_r = T_funR(Ti_vec_r);
sign_change_idx = find(T_fun_r(1:end-1).*T_fun_r(2:end) < 0);
Ti_zero_r = zeros(size(sign_change_idx));
for ki = 1:length(sign_change_idx)
    i = sign_change_idx(ki);
    % Valori ai due lati del cambio di segno
    x1 = Ti_vec_r(i);
    x2 = Ti_vec_r(i+1);
    y1 = T_fun_r(i);
    y2 = T_fun_r(i+1);
    % Linear interpolation
    Ti_zero_r(ki) = x1 - y1*(x2 - x1)/(y2 - y1);
end
%Choice of the supersonic solution
Ti_r = Ti_zero_r(1);
clear("Ti_zero_r","sign_change_idx","ki","x1",'"x2',"y1","y2")

% Rotor BC
Ti_r = fzero(T_funR,Ti_r);
Vi_r = sqrt(2*cp*(T0t-Ti_r));
Vit_r = V2t*(r2/ri_r);
Vin_r = sqrt(Vi_r^2 - Vit_r^2);
Mi_r = Vi_r/sqrt(gamma*R*Ti_r);
pi_r = p2*(Ti_r/T2)^(gamma/(gamma-1));
rhoi_r = pi_r/(R*Ti_r);
Wi_r = Vi_r - omega*ri_r;
Mwi_r = Wi_r/sqrt(gamma*R*Ti_r);


%% Plot

% Stator
figure
hold on
plot(Ti_vec_s,T_fun_s)
plot(Ti_s,T_fun(Ti_s),'ro')
yline(0,'--k')
xlabel('Stator Inlet Temperature')
ylabel('Tis - f(Tis)')
title('Stator Inlet BC')
%Blade Design Parameters
zeta_s = rad2deg(acos(0.5*(r0^2+r1^2-c_stator^2)/(r0*r1)));
Gamma_s = rad2deg(acos(0.5*(r1^2-r0^2+c_stator^2)/(c_stator*r1)));
psi_s = rad2deg(asin(sin(zeta_s)*r1/c_stator));

% Rotor
figure
hold on
plot(Ti_vec_r,T_fun_r)
plot(Ti_r,T_funR(Ti_r),'ro')
yline(0,'--k')
xlabel('Rotor Inlet Temperature')
ylabel('Tir - f(Tir)')
title('Rotor Inlet BC')



%% Stator and Rotor Hub/Shroud
if writeHubShroud == 1
    % HUB; SHROUD described in (r,theta,axial) (mm)
    HUB_s = [ri_s 0 0; ro_s 0 0]*1000;
    SHROUD_s = [ri_s 0 -b0;ro_s 0 -b0]*1000;
    % Write hub.crv
    fid = fopen('hub_s.crv','w');
    fprintf(fid, '%g %g %g\n', HUB_s.');   % 3 colonne
    fclose(fid);
    % Write shroud.crv
    fid = fopen('shroud_s.crv','w');
    fprintf(fid, '%g %g %g\n', SHROUD_s.');   % 3 colonne
    fclose(fid);



    HUB_r = [ri_r 0 0; r2 0 0; r3 0 0; ro_r 0 0]*1000;
    SHROUD_r = [ri_r 0 -b2; r2 0 -b2; r3 0 -b3; ro_r 0 -b3]*1000;
    % Write hub.crv
    fid = fopen('hub_r.crv','w');
    fprintf(fid, '%g %g %g\n', HUB_r.');   % 3 colonne
    fclose(fid);
    % Write shroud.crv
    fid = fopen('shroud_r.crv','w');
    fprintf(fid, '%g %g %g\n', SHROUD_r.');   % 3 colonne
    fclose(fid);

end