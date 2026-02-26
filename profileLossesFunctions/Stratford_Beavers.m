function [deltaS_sb,theta_sb,deltaSS_sb] = Stratford_Beavers(M0,M1,a1t,rho1t,mu1,c)
% Stratford Beavers Procedure
% 0: Blade Inlet
% 1: Blade Trailing Edge
% M: Mach Nummber
% at: Stagnation Sound Velocity
% rhot: Stagnation Density
% mu: Dynamic Viscosity
% c: Blade Chord
% errSB: if 1, it means Re is out of the range of the correlations. In this
% case BL properties are evaluated at Re = 10^6 and then Stewart losses
% will be corrected by a factor (Rex/Rex_ref)^-0.2 with Rex_ref = 10^6.

nu1t = mu1/rho1t;
P =  @(M) (M/(1+0.2*M^2))^4;
P1 = P(M1);
% Linear Function of M inside the channel
Pint =  @(x) ((M0+(M1-M0).*(x/c))./(1+0.2*(M0+(M1-M0).*(x/c)).^2)).^4;
X = integral(Pint,0,c)/P1;
% Rex evaluation
Rex = (a1t/nu1t)*X*M1*(1+0.2*M1^2)^(-2.5);

if (Rex > 5*10^5) && (Rex < 5*10^6)
    theta_sb = (0.036*(1+0.1*M1^2)^(-0.7))*X*Rex^(-1/5);
    deltaS_sb = 0.37*X*Rex^(-1/5);
    deltaSS_sb = (0.046*(1+0.8*M1^2)^0.44)*X*Rex^(-1/5);
elseif (Rex >= 5*10^6)
        theta_sb = (0.022*(1+0.1*M1^2)^(-0.7))*X*Rex^(-1/6);
        deltaS_sb = 0.23*X*Rex^(-1/6);
        deltaSS_sb = (0.028*(1+0.8*M1^2)^0.44)*X*Rex^(-1/6);
else
    % Re correction: all quantities are scaled by a factor (Re,ref/Re)^0.2
    K = (10^6/Rex)^0.2;
    theta_sb = K*(0.036*(1+0.1*M1^2)^(-0.7))*X*(10^6)^(-1/5);
    deltaS_sb = K*0.37*X*(10^6)^(-1/5);
    deltaSS_sb = K*(0.046*(1+0.8*M1^2)^0.44)*X*(10^6)^(-1/5);
end
end


