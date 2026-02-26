function [M2_stewart,ds_friction] = Stewart(t_te,q1,alpha1,M1_te,deltaS_sb,theta_sb,gamma,R,flagISO)
% t_te: trailing edge thickness
% q: blade pitch
% alpha: flow angle (relative if rotor) [rad]
% M: Mach number
% U: Peripheral Velocity
% deltaS_sb, thetaS_sb: boundary layer properties
% gamma, R: specific heta, gas constant
% flagISO: if 1 isentropic flow
% 0: Blade Inlet
% 1: Blade Trailing Edge
% 2: Blade Mixed Condition

% Since the procedure doesn't work so well close to meridional sonic
% conditions, I evaluate it for highly subsonic and highly supersonic flows
% and then interpolate the solution

if flagISO == 1
    M2_stewart = M1_te;
    ds_friction = 0;
else
    
    % Hp of the pitch increment due to radial configuration
    q2 = q1*1.05;
    % Evaluation of Stewart BL Thicknesses
    delta_te = t_te/(q1*cos(alpha1));
    deltaS = 2*deltaS_sb/(q1*cos(alpha1));
    theta = 2*theta_sb/(q1*cos(alpha1));
    M1_sonic_rad = 1/cos(alpha1);

    %% Do the procedure for the supersonic case changing K to have a resonable solution
    ds_friction = 1;
    i = 1;
    K = 1.5;
    while (ds_friction < 10) && (i < 500)
        M1 = K*M1_sonic_rad;
        %% Stewart Procedure
        V1_V1cr = M1*sqrt(0.5*(gamma+1)/(1+0.5*(gamma-1)*M1^2));
        A1fs = V1_V1cr^2*(gamma-1)/(gamma+1);
        rho1_rho1t = (1+0.5*(gamma-1)*M1^2)^(1/(1-gamma));

        C = ((1-A1fs)*0.5*(gamma+1)/gamma + cos(alpha1)^2*(1-deltaS-delta_te-theta)*V1_V1cr^2)/...
            (cos(alpha1)*(1-deltaS-delta_te)*V1_V1cr);

        D = V1_V1cr*sin(alpha1)*(1-deltaS-delta_te-theta)/(1-deltaS-delta_te);

        if M1*cos(alpha1) > 1
            % Supersonic Solution
            V2r_Vcr = gamma*C/(gamma+1) + sqrt((gamma*C/(gamma+1))^2 - 1 + D^2*(gamma-1)/(gamma+1));
        else
            % Subsonic Solution
            V2r_Vcr = gamma*C/(gamma+1) - sqrt((gamma*C/(gamma+1))^2 - 1 + D^2*(gamma-1)/(gamma+1));
        end
        rho2_rho2t = (1-(gamma-1)/(gamma+1)*(D^2 + V2r_Vcr^2))^(1/(gamma-1));
        p2t_p0t = (1-deltaS-delta_te)*cos(alpha1)*(V1_V1cr/V2r_Vcr)*(rho1_rho1t/rho2_rho2t)*(q1/q2);
        ds_friction = -R*log(p2t_p0t);

        K = K + 0.02;
        i = i + 1;
    end
    K = K - 0.02;

    M = [linspace(0.2,0.8*M1_sonic_rad,5) linspace(K*M1_sonic_rad,5,10)];
    ds_friction_plot = (1:length(M))*0;
    M2 = (1:length(M))*0;

    %% Atual Calculation
    for iM = 1:length(M)
        M1_plot = M(iM);
        %% Stewart Procedure
        V1_V1cr = M1_plot*sqrt(0.5*(gamma+1)/(1+0.5*(gamma-1)*M1_plot^2));
        A1fs = V1_V1cr^2*(gamma-1)/(gamma+1);
        rho1_rho1t = (1+0.5*(gamma-1)*M1_plot^2)^(1/(1-gamma));

        C = ((1-A1fs)*0.5*(gamma+1)/gamma + cos(alpha1)^2*(1-deltaS-delta_te-theta)*V1_V1cr^2)/...
            (cos(alpha1)*(1-deltaS-delta_te)*V1_V1cr);
        D = V1_V1cr*sin(alpha1)*(1-deltaS-delta_te-theta)/(1-deltaS-delta_te);

        if M1_plot*cos(alpha1) > 1
            % Supersonic Solution
            V2r_Vcr = gamma*C/(gamma+1) + sqrt((gamma*C/(gamma+1))^2 - 1 + D^2*(gamma-1)/(gamma+1));
        else
            % Subsonic Solution
            V2r_Vcr = gamma*C/(gamma+1) - sqrt((gamma*C/(gamma+1))^2 - 1 + D^2*(gamma-1)/(gamma+1));
        end
        rho2_rho2t = (1-(gamma-1)/(gamma+1)*(D^2 + V2r_Vcr^2))^(1/(gamma-1));
        p2t_p0t_plot = (1-deltaS-delta_te)*cos(alpha1)*(V1_V1cr/V2r_Vcr)*(rho1_rho1t/rho2_rho2t)*(q1/q2);
        M2r = V2r_Vcr*sqrt(2/(gamma+1-(gamma-1)*V2r_Vcr^2));
        M2(iM) = M2r/cos(alpha1);
        ds_friction_plot(iM) = -R*log(p2t_p0t_plot);
    end
    % Results Interpolation
   
    if M1_te < max(M)
        %ds_friction = spline(M,ds_friction_plot,M1_te);
        ds_friction = interp1(M,ds_friction_plot,M1_te, 'pchip');  %More stable solution
        M2_stewart = interp1(M,M2,M1_te, 'pchip');  %More stable solution
    else
        ds_friction = ds_friction_plot(end); 
        M2_stewart = M2(end); 
    end

end
end