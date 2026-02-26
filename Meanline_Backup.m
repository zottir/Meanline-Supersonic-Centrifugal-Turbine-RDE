clear;
close all;
clc;

flagPlot = 1; % if 1, plot
flagISO = 0;  % if 1, isentropic flow
flagMaps = 1; % if 1, not consider reflected shock

%% Data Definition
% Inlet Quantities Definition
T0 = 1500;
M0 = 2;
p0t = 10*10^5;
alpha0_deg = 0;
alpha0 = alpha0_deg*pi/180;
m = 1.5;

%Iterative Cycles
tolM3 = 10^-8;
NmaxM3 = 200;
tolTref = 10^-3;
NmaxTref = 100;
tolp1t = 10^-1;
Nmaxp1t = 100;
tolM1_stewart = 10^-8;
NmaxM1_stewart = 100;
tolp3t = 10^-1;
Nmaxp3t = 100;
tolMw3_stewart = 10^-8;
NmaxMw3_stewart = 100;

% Add to path all the subfolders of the folder shockFunctions
addpath(genpath('fluidPropertiesFunctions'));
addpath(genpath('profileLossesFunctions'));
addpath(genpath('shockLossesFunctions'));


% Fluid Properties Initializer
run("FluidPropInitializer.m");

phi = 2.5; %New phi definition as V1r/U2
lambda = 2;

r0 = 0.05; % Higashi used a cc of c.a 7 cm in diameter.
r1 = r0*1.5;
r2 = r1*1.05;
r3 = r2*1.32;
r4 = r3*1.05;

% Mean Radius
rm_stator = 0.5*(r0+r1);
rm_rotor = 0.5*(r2+r3);
% Meridional Chord
c_mer_stator = r1-r0;
c_mer_rotor = r3-r2;
% Chord
K_chord_s = 1.2;
K_chord_r = 1.1;
c_stator = K_chord_s*(r1-r0);
c_rotor = K_chord_r*(r3-r2);

% Incidence Angle [deg]
inc_LE_s = 0;
inc_LE_r = 0;
% Eccentricity LE
ecc_LE_s = 0;
ecc_LE_r = 0;

% Trailing Edges
t_TE_s = 0.5*10^-3; %0.5 mm
t_TE_r = 0.5*10^-3; %0.5 mm

% Zeifel Coefficients
C_Ft_s = 0.4;
C_Ft_r = 0.4;

% Peripheral Velocities
Umax = 500;
omega = Umax/r3;
U3 = Umax;
U2 = omega*r2;
rpm = omega/(2*pi/60);

% Flaring Angle
epsilonLim_deg = 5;
epsilonLim = epsilonLim_deg*pi/180;

%% Meanline Start
indTref = 0;
indp1t = 0;
indM1_te = 0;
indp3t = 0;
indM3 = 0;
indMw3_te = 0;


%First Guess for the reference temperature Tref = T0
Tref = T0;

for iTref = 1:NmaxTref
    % cp, gamma, R evaluation
    [gamma, R, MM, cp] = gamma_R_MM_cp_calc(Tref, dati, MolarFractions,MolecularMasses,zconst_low,zconst_up);
    % Viscosity evaluation
    [mu, k] = Sutherland(Tref,dati,dati_GasViscosityDatabase, ViscosityMixtureFormula);

    %% Stator Inlet
    % Calculation
    T0t = T0*(1+0.5*(gamma-1)*M0^2);
    rho0t = p0t/(R*T0t);
    p0 = p0t*(T0/T0t)^(gamma/(gamma-1));
    rho0 = p0/(R*T0);
    V0 = M0*sqrt(gamma*R*T0);
    V0r = V0*cos(alpha0);
    V0t = V0*sin(alpha0);
    a0 = sqrt(gamma*R*T0);
    a0t = sqrt(gamma*R*T0t);
    s0 = 0;

    b0 = m/(rho0*V0r*2*pi*r0);
    b1 = b0;
    b2 = b1;
    b3 = b2 - epsilonLim*(r3 - r2);
    b4 = b3 - epsilonLim*(r4 - r3);


    %% Stator Outlet
    T1t = T0t;
    m1 = m;
    V1r = phi*U2;
    rho1 = m1/(2*pi*r1*b1*V1r);

    % HP ISO: p1t = p0t
    p1t = p0t;
    for ip1t = 1:Nmaxp1t
        rho1t = p1t/(R*T1t);
        T1 = (rho1*R/p1t)^(gamma-1)*T1t^gamma;
        a1 = sqrt(gamma*R*T1);
        a1t = sqrt(gamma*R*T1t);
        p1 = rho1*R*T1;
        M1 = sqrt((T1t/T1 - 1)/(0.5*(gamma-1)));
        M1r = V1r/a1;
        V1 = M1*a1;
        V1t = sqrt(V1^2 - V1r^2);
        alpha1 = atan(V1t/V1r);
        alpha1_deg = alpha1*180/pi;
        DeflectionS = abs(alpha1_deg - alpha0_deg);

        % Geometry Definition
        % Number of Blades following Zweifel criteria for C_Ft_s = 0.4
        % Solidity defined along a reference mean radius.
        rad_solidity_stator = abs(2*cos(alpha1)^2*(V0r/V1r*tan(alpha0) - tan(alpha1))/C_Ft_s);
        q_m_stator = c_mer_stator/rad_solidity_stator;
        conc_s = 1;

        %Nb_stator = ceil(2*pi*rm_stator/q_m_stator);
        % Continuous to help convergence
        Nb_stator = 2*pi*rm_stator/q_m_stator;
        q1 = 2*pi*r1/Nb_stator;
        q0 = 2*pi*r0/Nb_stator;
        t_LE_s = q0/25;

        % Profile Losses Calculation
        % HP M1_te = M1
        M1_te = M1;
        for iM1_te = 1:NmaxM1_stewart
            [deltaS_sb_s,theta_sb_s,~] = Stratford_Beavers(M0,M1_te,a1t,rho1t,mu,c_stator);
            [M1_stewart,ds_friction_stator] = Stewart(t_TE_s,q1,alpha1,M1_te,deltaS_sb_s,theta_sb_s,gamma,R,flagISO);
            errM1_stewart(iM1_te + indM1_te(end)) = M1_stewart - M1;
            if abs(errM1_stewart(end)) < tolM1_stewart
                indM1_te = [indM1_te length(errM1_stewart)];
                break;
            else
                M1_te = M1_te*(1-0.4*errM1_stewart(end));
            end
        end

        % Shock Losses Calculation
        [shock_s,stators,phi_shock_s,inters_points_s,refl_s,ds_bowed_s,coll_shock_s] = bowShockLosses(M0,gamma,R,inc_LE_s,t_LE_s,ecc_LE_s,conc_s,alpha0,alpha1,Nb_stator,r0,r1,c_stator,flagISO);
        if flagMaps == 0
            [refl_shock_up_s,refl_shock_down_s,ds_refl_shock_s,detached_s,no_refl_Mn_s] = reflectShockLosses(M0,gamma,R,alpha0,r1,phi_shock_s,ds_bowed_s,stators,inters_points_s,refl_s,flagISO);
        else
            ds_refl_shock_s = 0;
        end
        clear inters_points_s;
        ds_shock_stator = ds_bowed_s + sum(ds_refl_shock_s);

        % Total Entropy Production in the Stator
        ds_stator = ds_friction_stator + ds_shock_stator;
        s1 = s0 + ds_stator;
        p1t_new = p0t*exp(-ds_stator/R);
        errp1t(ip1t + indp1t(end)) = p1t_new - p1t;
        if abs(errp1t(end)) < tolp1t
            indp1t = [indp1t length(errp1t)];
            break;
        else
            p1t = p1t_new;
        end
    end


    %% Rotor Inlet
    T2t = T1t;
    m2 = m1;
    V2t = V1t*(r1/r2);
    % ISO HP
    p2t = p1t;

    % Free-Vortex Mach function M = V/sqrt(gamma*R*T) forfsolve
    Mout_freeVortex = @(M, m, R, Tt, p_t, r, b, Vt, gamma) ...
        sqrt(((m*R*Tt./(2*pi*r*b*p_t)) .* (1+0.5*(gamma-1).*M.^2).^(1/(gamma-1))).^2 + Vt.^2) ...
        ./ sqrt(gamma*R*Tt ./ (1+0.5*(gamma-1).*M.^2));

    % Function to find M4
    fM2 = @(M2) M2 - Mout_freeVortex(M2, m2, R, T2t, p2t, r2, b2, V2t, gamma);
    % First Guess: Find the optimal zero depending on the supers/subs case
    xfindM2 = linspace(0.5*M1,2*M1,100);
    yfindM2 = fM2(xfindM2);
    indFindM2 = find(diff(sign(yfindM2)));
    if ~isempty(indFindM2) 
        if M1r > 1
            M2_guess = 0.5*( xfindM2(max(indFindM2)) + xfindM2(max(indFindM2)+1) );
        else
            M2_guess = 0.5*( xfindM2(min(indFindM2)) + xfindM2(min(indFindM2)+1) );
        end
    else
         M2_guess = M1;
    end
   
    options = optimset('TolX',1e-8);
    M2 = fzero(fM2, M2_guess, options);

    T2 = T2t/(1 + 0.5*(gamma-1)*M2^2);
    p2 = p2t*(T2/T2t)^(gamma/(gamma-1));
    rho2 = p2/(R*T2);
    a2 = sqrt(gamma*R*T2);
    a2t = sqrt(gamma*R*T2t);

    V2r = m/(2*pi*r2*b2*rho2);
    V2 = sqrt(V2r^2 + V2t^2);
    W2t = V2t - U2;
    W2r = V2r;
    W2 = sqrt(W2t^2 + W2r^2);
    M2r = V2r/a2;
    Mw2 = W2/a2;
    Mw2r = W2r/a2;
    alpha2 = atan(V2t/V2r);
    alpha2_deg = alpha2*180/pi;
    beta2 = atan(W2t/W2r);
    beta2_deg = beta2*180/pi;
    T2tw = T2 + 0.5*(W2^2-U2^2)/cp;
    p2tr = p2*(T2tw/T2)^(gamma/(gamma-1));
    s2 = s1 + cp*log(T2t/T1t) - R*log(p2t/p1t);


    %% Rotor Outlet
    m3 = m2;
    T3tw = T2tw;
    leu = lambda*U2^2;
    T3t = T2t - leu/cp;
    V3t = (V2t*U2 - leu)/U3;
    W3t = V3t - U3;

    %HP Iso
    p3t = p2t*(T3t/T2t)^(gamma/(gamma-1));
    for ip3t = 1:Nmaxp3t
        rho3t = p3t/(R*T3t);
        % Hp M3 = M2
        M3 = M2;
        for iM3 = 1:NmaxM3
            T3 = T3t/(1+0.5*(gamma-1)*M3^2);
            V3 = M3*sqrt(gamma*R*T3);
            V3r = sqrt(V3^2 - V3t^2);
            rho3 = m3/(2*pi*r3*b3*V3r);
            p3 = rho3*R*T3;
            M3_new = sqrt(((p3t/p3)^((gamma-1)/gamma) - 1)/(0.5*(gamma-1)));
            errM3(iM3 + indM3(end)) = M3_new - M3;
            if abs(errM3(end)) < tolM3
                indM3 = [indM3 length(errM3)];
                break;
            else
                M3 = 0.7*M3 +0.3*M3_new;
            end
        end

        % fun = @(M3) M3 - M3_new_from_M3(M3, gamma, R, T3t, p3t, V3t, m3, r3, b3);
        % % First Guess
        % M3_guess = M2;  
        % M3 = fzero(fun, M3_guess, options);
        % 
        % T3 = T3t/(1+0.5*(gamma-1)*M3^2);
        % V3 = M3*sqrt(gamma*R*T3);
        % V3r = sqrt(V3^2 - V3t^2);
        % rho3 = m3/(2*pi*r3*b3*V3r);
        % p3 = rho3*R*T3;

        W3r = V3r;
        W3 = sqrt(W3t^2 + W3r^2);
        a3 = sqrt(gamma*R*T3);
        a3t = sqrt(gamma*R*T3t);
        Mw3 = W3/a3;
        Mw3r = W3r/a3;
        M3r = V3r/a3;
        p3tr = p3*(T3tw/T3)^(gamma/(gamma-1));
        alpha3 = atan(V3t/V3r);
        alpha3_deg = alpha3*180/pi;
        beta3 = atan(W3t/W3r);
        beta3_deg = beta3*180/pi;
        DeflectionR = abs(beta3_deg - beta2_deg);

        % Geometry Definition
        % Number of Blades following Zweifel criteria for C_Ft = 0.4
        % Solidity defined along a reference mean radius.
        rad_solidity_rotor = abs(2*cos(beta2)^2*(V2r/V3r*tan(beta2) - tan(beta3))/C_Ft_r);
        q_m_rotor = c_mer_rotor/rad_solidity_rotor;
        conc_r = 1;
        %Nb_rotor = ceil(2*pi*rm_rotor/q_m_rotor);
        Nb_rotor = 2*pi*rm_rotor/q_m_rotor; % Continuous to help convergence

        q2 = 2*pi*r2/Nb_rotor;
        q3 = 2*pi*r3/Nb_rotor;
        t_LE_r = q2/25;

        % Profile Losses Calculation
        % HP Mw3_te = Mw3
        Mw3_te = Mw3;
        for iMw3_te = 1:NmaxMw3_stewart
            [deltaS_sb_r,theta_sb_r,~] = Stratford_Beavers(Mw2,Mw3_te,a3t,rho3t,mu,c_rotor);
            [Mw3_stewart,ds_friction_rotor] = Stewart(t_TE_r,q3,beta3,Mw3_te,deltaS_sb_r,theta_sb_r,gamma,R,flagISO);
            errMw3_stewart(iMw3_te + indMw3_te(end)) = Mw3_stewart - Mw3;
            if abs(errMw3_stewart(end)) < tolMw3_stewart
                indMw3_te = [indMw3_te length(errMw3_stewart)];
                break;
            else
                Mw3_te = Mw3_te * (1 - 0.4*errMw3_stewart(end)/Mw3);
            end
        end
       

        % Shock Losses Calculation
        if Mw2 > 1.15
            [shock_r,rotors,phi_shock_r,inters_points_r,refl_r,ds_bowed_r,coll_shock_r] = bowShockLosses(Mw2,gamma,R,inc_LE_r,t_LE_r,ecc_LE_r,conc_r,beta2,beta3,Nb_rotor,r2,r3,c_rotor,flagISO);
            if flagMaps == 0
                [refl_shock_up_r,refl_shock_down_r,ds_refl_shock_r,detached_r,no_refl_Mn_r] = reflectShockLosses(Mw2,gamma,R,beta2,r3,phi_shock_r,ds_bowed_r,rotors,inters_points_r,refl_r,flagISO);
            else
                ds_refl_shock_r = 0;
            end
            clear inters_points_r;
            ds_shock_rotor = ds_bowed_r + sum(ds_refl_shock_r);
        else
            [rotors] = approximate_blade(beta2,beta3,r2,r3,c_rotor,Nb_rotor,conc_r);
            ds_shock_rotor = 0;
            refl_shock_up_r = [];
            refl_shock_down_r = [];
            shock_r = [];
            coll_shock_r = 0;
            detached_r = 0;
            no_refl_Mn_r = 0;
        end

        ds_rotor = ds_friction_rotor + ds_shock_rotor;
        s3 = s2 + ds_rotor;
        p3t_new = p2t*exp(cp/R*log(T3t/T2t) - ds_rotor/R);
        errp3t(ip3t + indp3t(end)) = p3t_new - p3t;
        if abs(errp3t(end)) < tolp3t
            indp3t = [indp3t length(errp3t)];
            break;
        else
            p3t = p3t_new;
        end
    end

    % Reference Temperature Update for the fluid properties
    Tref_new = 0.5*(max([T0,T1,T2,T3]) + min([T0,T1,T2,T3]));
    errTref(iTref + indTref(end)) = Tref_new - Tref;
    if abs(errTref(end)) < tolTref
        indTref = [indTref length(errTref)];
        break;
    end
    Tref = Tref_new;
end


%% Machine Outlet
T4t = T3t;
m4 = m3;
V4t = V3t*(r3/r4);
% ISO HP
p4t = p3t;

% Function to find M4
fM4 = @(M4) M4 - Mout_freeVortex(M4, m4, R, T4t, p4t, r4, b4, V4t, gamma);
% First Guess: Find the optimal zero depending on the supers/subs case
    xfindM4 = linspace(0.5*M3,2*M3,100);
    yfindM4 = fM4(xfindM4);
    indFindM4 = find(diff(sign(yfindM4)));
    if ~isempty(indFindM4) 
        if M3r > 1
            M4_guess = 0.5*( xfindM4(max(indFindM4)) + xfindM4(max(indFindM4)+1) );
        else
            M4_guess = 0.5*( xfindM4(min(indFindM4)) + xfindM4(min(indFindM4)+1) );
        end
    else
         M4_guess = M3;
    end
M4 = fzero(fM4, M4_guess, options);

T4 = T4t / (1 + 0.5*(gamma-1)*M4^2);
p4 = p4t * (T4/T4t)^(gamma/(gamma-1));
rho4 = p4 / (R*T4);
V4r = m / (2*pi*r4*b4*rho4);
V4 = sqrt(V4r^2 + V4t^2);
a4 = sqrt(gamma * R * T4);
a4t = sqrt(gamma*R*T4t);
M4r = V4r/a4;
alpha4 = atan(V4t/V4r);
alpha4_deg = alpha4*180/pi;
s4 = s3 + cp*log(T4t/T3t) - R*log(p4t/p3t);


T4t_is = T0t*(p4t/p0t)^((gamma-1)/gamma);
efficiency_tt = (T0t - T4t)/(T0t - T4t_is);
Power = m*leu;


%% Check
% Mass Flows
m1_test = 2*pi*r1*b1*rho1*V1r;
m2_test = 2*pi*r2*b2*rho2*V2r;
m3_test = 2*pi*r3*b3*rho3*V3r;
if (abs(m1_test-m1) > 10^-6)
    fprintf('Mass Flow Rate Conservation NOT RESPECTED at Stator Outlet:\n m0 = %d \n m1 = %d \n',m,m1);
end
if (abs(m2_test-m2) > 10^-6)
    fprintf('Mass Flow Rate Conservation NOT RESPECTED at Rotor Inlet:\n m1 = %d \n m2 = %d \n',m1,m2);
end
if (abs(m3_test-m3) > 10^-6)
    fprintf('Mass Flow Rate Conservation NOT RESPECTED at Rotor Outlet:\n m2 = %d \n m3 = %d \n',m2,m3);
end


% Rothalpy Check
Rotalphy2 = cp*T2 +0.5*(W2^2 - U2^2);
Rotalphy3 = cp*T3 +0.5*(W3^2 - U3^2);
if (abs(Rotalphy2-Rotalphy3) > 10^-1)
    fprintf('Rothalphy Conservation NOT RESPECTED\n');
end


% Check 1 < M1 < M2
if (M1r<1 && M2r>M1r )
    fprintf('Unphysical solution: M1r = %.2f < 1 cant lead to M2r = %.2f > M1r\n',M1r,M2r);
end
if (M1r>1 && M2r<M1r )
    fprintf('Unphysical solution: M1r = %.2f > 1 cant lead to M2r = %.2f < M1r\n',M1r,M2r);
end

% Check Kantrowitz
M_plot = linspace(1,5,1000);
Mk_plot = sqrt(((gamma-1)*M_plot.^2+2)./(2*gamma*M_plot.^2-(gamma-1)));
ARis_plot = M_plot.*(0.5*(gamma+1)./(1+0.5*(gamma-1)*M_plot.^2)).^(0.5*(gamma+1)/(gamma-1));
AR_plot = Mk_plot.*(0.5*(gamma+1)./(1+0.5*(gamma-1)*Mk_plot.^2)).^(0.5*(gamma+1)/(gamma-1));
% Stator
Kantrowitz_s = 0;
Mk_s = sqrt(((gamma-1)*M0^2+2)/(2*gamma*M0^2-(gamma-1)));
AR_min_s = Mk_s*(0.5*(gamma+1)/(1+0.5*(gamma-1)*Mk_s^2))^(0.5*(gamma+1)/(gamma-1));
AR_s = (b1/b0)*(r1/r0)*(cos(alpha1)/cos(alpha0));
if (AR_s < AR_min_s)
    Kantrowitz_s = 1;
    fprintf('Kantowitz limit NOT respected for stator\n');
end
if coll_shock_s == 1
    fprintf('Collective Shock in the Stator\n');
end
% Rotor
Kantrowitz_r = 0;
AR_r = (b3/b2)*(r3/r2)*(cos(beta3)/cos(beta2));
Mk_r = sqrt(((gamma-1)*Mw2^2+2)/(2*gamma*Mw2^2-(gamma-1)));
AR_min_r = Mk_r*(0.5*(gamma+1)/(1+0.5*(gamma-1)*Mk_r^2))^(0.5*(gamma+1)/(gamma-1));
if (AR_r < AR_min_r) && (Mw2 >1)
    Kantrowitz_r = 1;
    fprintf('Kantowitz limit NOT respected for rotor\n');
end
if coll_shock_r == 1
    fprintf('Collective Shock in the Rotor\n');
end


% % Info Rotor Shape
% if (Mw2 < 1) && (Mw3 < Mw2)
%     fprintf('Mw2 = %.2f and Mw3 = %.2f : Divergent rotor shape to decelerate the flow\n',Mw2,Mw3);
% elseif (Mw2 < 1) && (Mw3 > 1)
%     fprintf('Mw2 = %.2f and Mw3 = %.2f : Convergent-Divergent rotor shape to accelerate the flow\n',Mw2,Mw3);
% elseif (Mw2 < 1) && (Mw3 > Mw2) && (Mw3 < 1)
%     fprintf('Mw2 = %.2f and Mw3 = %.2f : Convergent rotor shape to accelerate the flow\n',Mw2,Mw3);
% elseif (Mw2 > 1) && (Mw3 > Mw2)
%     fprintf('Mw2 = %.2f and Mw3 = %.2f : Divergent rotor shape to accelerate the flow\n',Mw2,Mw3);
% elseif (Mw2 > 1) && (Mw3 < 1)
%     fprintf('Mw2 = %.2f and Mw3 = %.2f : Convergent-Divergent rotor shape to decelerate the flow\n',Mw2,Mw3);
% elseif (Mw2 > 1) && (Mw3 < Mw2) && (Mw3 > 1)
%     fprintf('Mw2 = %.2f and Mw3 = %.2f : Convergent rotor shape to decelerate the flow\n',Mw2,Mw3);
% end


% % Info Reflected Shock
% if flagMaps == 0
%     if detached_s(1) == 1
%         fprintf('Stator''s upper reflected shock is detached\n');
%     end
%     if detached_s(2) == 1
%         fprintf('Stator''s lower reflected shock is detached\n');
%     end
%     if detached_r(1) == 1
%         fprintf('Rotor''s upper reflected shock is detached\n');
%     end
%     if detached_r(2) == 1
%         fprintf('Rotor''s lower reflected shock is detached\n');
%     end
% 
%     if no_refl_Mn_s(1) == 1
%         fprintf('Stator''s upper reflected shock not present (M after bow shock < 1)\n');
%     end
%     if no_refl_Mn_s(2) == 1
%         fprintf('Stator''s lower reflected shock not present (M after bow shock < 1)\n');
%     end
%     if no_refl_Mn_r(1) == 1
%         fprintf('Rotor''s upper reflected shock not present (M after bow shock < 1)\n');
%     end
%     if no_refl_Mn_r(2) == 1
%         fprintf('Rotor''s lower reflected shock not present (M after bow shock < 1)\n');
%     end
% end


%% Plot
if flagPlot == 1

    % Latex Format
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');


    % T - s Diagram
    figure
    hold on;
    plot(s0,T0t,'o',s1,T1t,'o',s2,T2t,'o',s3,T3t,'o',s4,T4t,'o', 'MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(s0,T0,'o',s1,T1,'o',s2,T2,'o',s3,T3,'o',s4,T4,'o', 'MarkerFaceColor','b','MarkerEdgeColor','b');
    plot(s2,T2tw,'o',s3,T3tw,'o', 'MarkerFaceColor','r','MarkerEdgeColor','r');
    plot(0-max([s0 s1 s2 s3 s4])/5, min([T0,T1,T2,T3,T4]) - (T0t - min([T0,T1,T2,T3,T4]))/4);
    plot(max([s0 s1 s2 s3 s4])+max([s0 s1 s2 s3 s4])/4,T0t + (T0t - min([T0,T1,T2,T3,T4]))/4);

    plot([s0 s1 s2 s3],[T0t T1t T2t T3t],'k','LineWidth',1.2);
    plot([s0 s1 s2 s3],[T0 T1 T2 T3],'b','LineWidth',1.2);
    plot([s2 s3],[T2tw T3tw],'r','LineWidth',1.2);

    textOpts = {'FontSize',14,'FontWeight','bold','VerticalAlignment','bottom','HorizontalAlignment','left'};
    text(s0,T0,'$T_0$',textOpts{:});
    text(s1,T1,'$T_1$',textOpts{:});
    text(s2,T2,'$T_2$',textOpts{:});
    text(s3,T3,'$T_3$',textOpts{:});
    text(s4,T4,'$T_4$',textOpts{:});
    text(s0,T0t,'$T_{0T}$',textOpts{:});
    text(s1,T1t,'$T_{1T} = T_{2T}$',textOpts{:});
    text(s3,T3t,'$T_{3T} = T_{4T}$',textOpts{:});
    text(s2,T2tw,'$T_{2Tr}$',textOpts{:});
    text(s3,T3tw,'$T_{3Tr}$',textOpts{:});

    xlabel('Specific Entropy $s$ [J/kgK]');
    ylabel('Temperature T [K]');
    set(gca,'FontSize',14);
    grid on;


    % Kantrowitz Limit
    figure
    hold on;
    plot(M_plot,ARis_plot,'k',M_plot,AR_plot,'b','LineWidth',1.2);
    plot(M0,AR_s,'ro','MarkerFaceColor','r','MarkerEdgeColor','r');
    plot(Mw2,AR_r,'ro','MarkerFaceColor','r','MarkerEdgeColor','r');
    text(M0 + 0.1,AR_s,'S','FontSize',14,'FontWeight','bold');
    text(Mw2 + 0.1,AR_r,'R','FontSize',14,'FontWeight','bold');
    xlabel('$M_{in}$');
    ylabel('$A_{out}/A_{in}$');
    title('Kantrowitz Limit');
    set(gca,'FontSize',14);
    grid on;


    % Velocity Triangles
    % quiver(x_coda,y_coda, V(1), V(2), no autoscale)
    arancio = [1 0.8 0];
    arancio_scuro = [1 0.5 0];
    rosso = [1 0 0];
    viola = [0.5 0 0.5];
    figure
    %Stator
    subplot(1,3,1);
    hold on;
    quiver(0, 0, V0t, V0r, 0,'Color',viola,'LineWidth',1.2) % V0
    quiver(0, 0, V1t, V1r, 0,'Color',rosso,'LineWidth',1.2) % V1
    xlim([min([0,V0t,V1t,U2,V2t,W2t,U3,V3t,W3t])-100,max([0,V0t,V1t,U2,V2t,W2t,U3,V3t,W3t])+50]);
    ylim([-100,max([V0r,V1r,V2r,V3r])+50]);
    legend('$V_0$','$V_1$');
    ylabel('Radial Component [m/s]');
    title('Stator');
    grid on;
    set(gca,'FontSize',14);
    % Rotor Inlet
    subplot(1,3,2);
    hold on;
    quiver(0, 0, U2, 0, 0, 'Color', arancio,'LineWidth',1.2) % U2
    quiver(U2, 0, W2t, V2r, 0,'Color', arancio_scuro,'LineWidth',1.2) % W2
    quiver(0, 0, V2t, V2r, 0,'Color',rosso,'LineWidth',1.2) % V2
    xlim([min([0,V0t,V1t,U2,V2t,W2t,U3,V3t,W3t])-100,max([0,V0t,V1t,U2,V2t,W2t,U3,V3t,W3t])+50]);
    ylim([-100,max([V0r,V1r,V2r,V3r])+50]);
    legend('$U_2$','$W_2$','$V_2$');
    xlabel('Tangential Component [m/s]');
    title('Rotor Inlet');
    grid on;
    set(gca,'FontSize',14);
    % Rotor Outlet
    subplot(1,3,3)
    hold on;
    quiver(0, 0, U3, 0, 0,'Color', arancio,'LineWidth',1.2) % U3
    quiver(U3, 0, W3t, V3r, 0,'Color',arancio_scuro,'LineWidth',1.2) % W3
    quiver(0, 0, V3t, V3r, 0,'Color',rosso,'LineWidth',1.2) % V3
    xlim([min([0,V0t,V1t,U2,V2t,W2t,U3,V3t,W3t])-100,max([0,V0t,V1t,U2,V2t,W2t,U3,V3t,W3t])+50]);
    ylim([-100,max([V0r,V1r,V2r,V3r])+50]);
    legend('$U_3$','$W_3$','$V_3$');
    title('Rotor Outlet');
    grid on;
    set(gca,'FontSize',14);
    % % Rotor Velocity Triangle
    % figure
    % hold on;
    % quiver(0, 0, U2, 0, 0,'LineWidth',1.2) % U2
    % quiver(U2, 0, W2t, V2r, 0,'LineWidth',1.2) % W2
    % quiver(0, 0, V2t, V2r, 0,'LineWidth',1.2) % V2
    % quiver(0, 0, U3, 0, 0,'LineWidth',1.2) % U3
    % quiver(U3, 0, W3t, V3r, 0,'LineWidth',1.2) % W3
    % quiver(0, 0, V3t, V3r, 0,'LineWidth',1.2) % V3
    % xlim([min([0,V0t,V1t,U2,V2t,W2t,U3,V3t,W3t])-100,max([0,V0t,V1t,U2,V2t,W2t,U3,V3t,W3t])+50]);
    % ylim([-100,max([V0r,V1r,V2r,V3r])+50]);
    % legend('$U_2$','$W_2$','$V_2$','$U_3$','$W_3$','$V_3$');
    % xlabel('Tangential Component [m/s]');
    % ylabel('Radial Component [m/s]');
    % title('Rotor');
    % grid on;
    % set(gca,'FontSize',14);

    % % Solution fsolve free-vortex
    % M2_range = linspace(0.1, 3, 200);
    % f_vals = arrayfun(fM2, M2_range);
    % figure;
    % plot(M2_range, f_vals, 'LineWidth', 1.5);
    % hold on; grid on;
    % yline(0,'--k');
    % xline(M2,'--r','M2 solution');
    % xline(M1,'--r','M1');
    % xlabel('$M_2$');
    % ylabel('$f(M_2) = M_2 - V_2/a_2$');
    % title('Free Vortex Region M Research');
    % set(gca,'FontSize',14);

    % M4_range = linspace(0.1, 3, 200);
    % f_vals = arrayfun(fM4, M4_range);
    % figure;
    % plot(M4_range, f_vals, 'LineWidth', 1.5);
    % hold on; grid on;
    % yline(0,'--k');
    % xline(M4,'--r','M4 ');
    % xline(M3,'--r','M3');
    % xlabel('$M_3$');
    % ylabel('$f(M_4) = M_4 - V_4/a_4$');
    % title('Free Vortex Region M Research');
    % set(gca,'FontSize',14);


    % Shock Patterns
    figure
    hold on
    % Radii
    circ = linspace(0,2*pi,500);
    plot(r0*cos(circ),r0*sin(circ),'g')
    plot(r1*cos(circ),r1*sin(circ),'g')
    plot(r2*cos(circ),r2*sin(circ),'g')
    plot(r3*cos(circ),r3*sin(circ),'g')
    % Shocks
    plot(shock_s(1,:),shock_s(2,:),'r','LineWidth',1);
    if ~isempty(shock_r )
    plot(shock_r(1,:),shock_r(2,:),'r','LineWidth',1);
    end
    if flagMaps == 0
        if ~isempty(refl_shock_up_s)
            plot(refl_shock_up_s(1,:),refl_shock_up_s(2,:),'r','LineWidth',1);
        end
        if ~isempty(refl_shock_down_s)
            plot(refl_shock_down_s(1,:),refl_shock_down_s(2,:),'r','LineWidth',1);
        end
        if ~isempty(refl_shock_up_r)
            plot(refl_shock_up_r(1,:),refl_shock_up_r(2,:),'r','LineWidth',1);
        end
        if ~isempty(refl_shock_down_r)
            plot(refl_shock_down_r(1,:),refl_shock_down_r(2,:),'r','LineWidth',1);
        end
    end
    % Approximate blades
    plot(stators(1,:),stators(2,:),'k',stators(3,:),stators(4,:),'k',stators(5,:),stators(6,:),'k','LineWidth',1.25);
    plot(rotors(1,:),rotors(2,:),'k',rotors(3,:),rotors(4,:),'k',rotors(5,:),rotors(6,:),'k','LineWidth',1.25);
    % Settings
    grid on;
    axis equal;
    set(gca,'FontSize',14);
    xlabel('x [m]');
    ylabel('y [m]');


    % Convergence
    figure
    plot(abs(errTref))
    ylabel('$error T_{ref}$')
    xlabel('Iterations [-]')
    grid on
    set(gca,'FontSize',14);
    set(gca,'YScale','log');

    figure
    hold on
    plot(abs(errp1t))
    xline(indp1t,'--r');
    ylabel('$error p_{1t}$')
    xlabel('Iterations [-]')
    grid on
    set(gca,'FontSize',14);
    set(gca,'YScale','log');

    figure
    hold on
    plot(abs(errp3t));
    xline(indp3t,'--r');
    ylabel('$error p_{3t}$')
    xlabel('Iterations [-]')
    grid on
    set(gca,'FontSize',14);
    set(gca,'YScale','log');

    figure
    hold on
    plot(abs(errM3))
    xline(indM3,'--r');
    ylabel('$error M_3$')
    xlabel('Iterations [-]')
    grid on
    set(gca,'FontSize',14);
    set(gca,'YScale','log');

    figure
    hold on;
    plot(abs(errM1_stewart))
    xline(indM1_te,'--r');
    ylabel('$error M_{1,Stewart}$')
    xlabel('Iterations [-]')
    grid on
    set(gca,'FontSize',14);
    set(gca,'YScale','log');

    figure
    hold on
    plot(abs(errMw3_stewart))
    xline(indMw3_te,'--r');
    ylabel('$error M_{w3}$')
    xlabel('Iterations [-]')
    grid on
    set(gca,'FontSize',14);
    set(gca,'YScale','log');

end




