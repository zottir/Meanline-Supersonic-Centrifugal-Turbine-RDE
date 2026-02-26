function [deltax,x_shock_up,y_shock_up,x_shock_down, y_shock_down,...
    phi_up,phi_down,x0,SB,S] = bowed_shockWave...
    (M_in,t_LE,x_stat,y_stat,alpha_g, i_angle, ecc_LE, gamma, R, phiS_ratio, geom, optionFinalSlopeCorrect)
% Function to draw the bowed shock ahead of the leading edge
% INPUT:
% M_in: Inlet Mach Number;
% t_LE: Leading Edge Thickness [m];
% x_stat, y_stat: Coordinated of the rounded elliptical leading edge;
% alpha_g: Leading-edge geometric angle [rad];
% i_angle: Leading-edge incidence angle  [rad];
% ecc_le: Eccentricity of the leading Edge [-];
% gamma: Specific heat ratio [-];
% R: Gas Constant;          
% phiS_ratio: Geometrical parameter equal to 6;
% geom: 1 if 2D geometry, 2 if axial-symmetric
% optionFinalSlopeCorrect: if 1 correct the final slope of the bowed shock
%
% OUTPUT
% deltax: Distance between Shock foremost point and LE;
% x_shock, y_shock: bow shock coordinates
% phi_up, phi_down: shock angle [rad]
% x0: Shock Foremost point
% SB: Body Sonic Point Coordinates
% |x_sb_down x_sb_up |
% |y_sb_down y_sb_up |
% S: Sonic Point Coordinates
% |x_s_down x_s_up |
% |y_s_down y_s_up |

%% Preliminary Calculus
%ridefinisco alpha per tenere conto dell'angolo di incidenza
alpha = alpha_g + i_angle;

%% Moeckel
[~, phi_C, Ps_by_P0_c, C_shock, sigma_shock, lambda_S] = Moeckel(M_in, gamma, R, geom);
% Definition of thetad by points
M_nodi = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 3.0 3.4];
theta_d_nodi = [0 21 28 34 38 41 43 45.5 46.5 49 50.5];
theta_d = deg2rad(cubicspline(M_nodi, theta_d_nodi, M_in));
phi_S=phi_C;
lamda_d = lambda_S;

beta_shock = sqrt(M_in^2 - 1);
B = sigma_shock/Ps_by_P0_c;

[beta_phiS_down] = betaBowShockModel(M_in, phiS_ratio, gamma, R);
if (optionFinalSlopeCorrect == 1) && isreal(beta_phiS_down)
    %Calcolo pendenza finale del bow shock utilizzando il modello che si basa
    %sull'analogia con un urto obliquo da un corner
    m_shock_inf_down = tan(beta_phiS_down);
    beta_shock = 1/m_shock_inf_down;
end

if geom == 1      % geometria 2D
    eta = (lambda_S + lamda_d)/2;
    yS_ySB = (1 - B*cos(eta))^(-1);
elseif geom == 2    % geometria assial simmetrica
    eta = (lambda_S + theta_d)/2;
    yS_ySB = (1 - B*cos(eta))^(-1/2);
end

%% Profilo leading edge ellisse
R_sh = (t_LE)/2;
M = [0; R_sh];
N = [0; -R_sh];

C_LE = [(M(1,1)+N(1,1))/2; (M(2,1)+N(2,1))/2];
b_LE = (M(2,1)-N(2,1))/2;
a_LE = sqrt((b_LE^2)/(1-ecc_LE^2));

x_LE = linspace(C_LE(1,1)-a_LE, C_LE(1,1),1500);
%upper part of l.e.
y_LE_up = C_LE(2,1)+sqrt((1 - (x_LE-C_LE(1,1)).^2./(a_LE.^2))*(b_LE^2));
x_LE_up = x_LE(end:-1:1);
y_LE_up = y_LE_up(end:-1:1);
%lower part of l.e. 
y_LE_down=C_LE(2,1)-sqrt((1 - (x_LE-C_LE(1,1)).^2./(a_LE.^2))*(b_LE^2));
x_LE_down = x_LE;

%applico la rotazione per tenere conto dell'angolo di incidenza  
[x_LE_up, y_LE_up] = rotation(x_LE_up, y_LE_up, -i_angle);
[x_LE_down, y_LE_down] = rotation(x_LE_down, y_LE_down, -i_angle);

%Calcolo pendenza locale per trovare il sonic point on body 
m_LE_up = (y_LE_up(2:end) - y_LE_up(1:end-1))./(x_LE_up(2:end) - x_LE_up(1:end-1));
m_LE_down = (y_LE_down(2:end) - y_LE_down(1:end-1))./(x_LE_down(2:end) - x_LE_down(1:end-1));
[~,I_SB_up] = min(abs(atan(m_LE_up)-eta));
[~,I_SB_down] = min(abs((atan(m_LE_down))+eta));
x_SB_up = x_LE_up(I_SB_up);
y_SB_up = y_LE_up(I_SB_up);
x_SB_down = x_LE_down(I_SB_down);
y_SB_down = y_LE_down(I_SB_down);

%Troviamo il punto di intersezione delle tangenti nei SB
m_up = tan(eta);
y_tan0 = (m_up*(x_SB_down - x_SB_up) + (y_SB_up + y_SB_down))/(2);

%Shift tutte le curve per porre la streamline di ristagno ad altezza nulla
y_SB_up = y_SB_up - y_tan0;
y_SB_down = y_SB_down - y_tan0;

%% Forma bow shock upper 

x_0_up = y_SB_up*beta_shock*yS_ySB*sqrt(((beta_shock*tan(phi_S))^2)-1);
L_up = y_SB_up*((yS_ySB*(C_shock + tan(eta))) - tan(eta));
x_SB_up = x_0_up + L_up;

%traslo ellisse per porlo nel sistema di riferimento dello schock
delta_up = (x_LE_up(I_SB_up)-x_SB_up);
x_LE_up = x_LE_up - delta_up;

%punto sonico
x_S_up = (beta_shock*x_0_up)/(sqrt((beta_shock)^2 - ((tan(phi_S))^(-2))));
y_S_up = (x_0_up)/(beta_shock*tan(phi_S)*sqrt((beta_shock^2) - ((tan(phi_S))^(-2))));

%shock
x_g_up = sqrt((beta_shock*y_S_up)^2 + x_0_up^2);
x_shock_up = linspace(x_0_up, 10*x_g_up, 500);
y_shock_up = (sqrt(((x_shock_up).^2) - ((x_0_up).^2)))./(beta_shock);
dydx_up = (sqrt(x_0_up.^2 + ((beta_shock).^2.*(y_shock_up).^2)))./(((beta_shock).^2).*y_shock_up);


[beta_phiS_up] = betaBowShockModel(M_in, phiS_ratio, gamma, R);
if (optionFinalSlopeCorrect == 1) && isreal(beta_phiS_up)
    %Calcolo pendenza finale del bow shock utilizzando il modello che si basa
    %sull'analogia con un urto obliquo da un corner
    m_shock_inf_up = tan(beta_phiS_up);

    %Individuo la pendenza corrispondente a m_shock_inf
    [~, I_minf_up] = min(abs(dydx_up - m_shock_inf_up));
    m_shock_up = dydx_up(I_minf_up);
    x_shock_inf_up = linspace(x_shock_up(I_minf_up), 12*x_g_up, 500);
    y_shock_inf_up = y_shock_up(I_minf_up) + m_shock_up*(x_shock_inf_up - x_shock_up(I_minf_up));

    x_shock_up = [x_shock_up(1:I_minf_up) x_shock_inf_up];
    y_shock_up = [y_shock_up(1:I_minf_up) y_shock_inf_up];

end

%% Forma bow shock lower 

%Rifletto la parte lower per utilizzare le stesse formule e non cambiare i
%segni
y_SB_down = -y_SB_down;

x0_down = y_SB_down*beta_shock*yS_ySB*sqrt(((beta_shock*tan(phi_S))^2)-1);
L_down = y_SB_down*((yS_ySB*(C_shock + tan(eta))) - tan(eta));
x_SB_down = x0_down + L_down;

%traslo ellisse per porlo nel sistema di riferimento dello schock
delta_down = (x_LE_down(I_SB_down)-x_SB_down);
x_LE_down = x_LE_down - delta_down;

%punto sonico
x_S_down = (beta_shock*x0_down)/(sqrt((beta_shock)^2 - ((tan(phi_S))^(-2))));
y_S_down = (x0_down)/(beta_shock*tan(phi_S)*sqrt((beta_shock^2) - ((tan(phi_S))^(-2))));

%shock
x_g_down = sqrt((beta_shock*y_S_down)^2 + x0_down^2);
x_shock_down = linspace(x0_down, 10*x_g_down, 500);
y_shock_down = (sqrt(((x_shock_down).^2) - ((x0_down).^2)))./(beta_shock);
dydx_down = (sqrt(x0_down.^2 + ((beta_shock).^2.*(y_shock_down).^2)))./(((beta_shock).^2).*y_shock_down);

[beta_phiS_down] = betaBowShockModel(M_in, phiS_ratio, gamma, R);
if (optionFinalSlopeCorrect == 1) && isreal(beta_phiS_down)
    %Calcolo pendenza finale del bow shock utilizzando il modello che si basa
    %sull'analogia con un urto obliquo da un corner
    m_shock_inf_down = tan(beta_phiS_down);

    %Individuo la pendenza corrispondente a m_shock_inf
    [~, I_minf_down] = min(abs(dydx_down - m_shock_inf_down));
    m_shock_down = dydx_down(I_minf_down);
    x_shock_inf_down = linspace(x_shock_down(I_minf_down), 12*x_g_down, 500);
    y_shock_inf_down = y_shock_down(I_minf_down) + m_shock_down*(x_shock_inf_down - x_shock_down(I_minf_down));
    x_shock_down = [x_shock_down(1:I_minf_down) x_shock_inf_down];
    y_shock_down = [y_shock_down(1:I_minf_down) y_shock_inf_down];

end

%% Rifletto per tornare al sistema di riferimento originale
y_S_down = -y_S_down;
y_shock_down = -y_shock_down;
y_SB_down = -y_SB_down;

%% Shift per continuità del profilo della l.e.
deltax_shift = x_LE_down(1) - x_LE_up(end);
x_LE_down = x_LE_down - deltax_shift;
x_SB_down = x_SB_down - deltax_shift;
x_S_down = x_S_down - deltax_shift;
x_shock_down = x_shock_down - deltax_shift;

%% Rotazione in funzione dell'angolo di ingresso
[x_stat,~] = rotation(x_stat', y_stat', -alpha);   %angolo negativo perché devo ruotare i punti in senso orario
%shift dei profili per far coincidere le leading edge (individuo i due punti corrispondenti)
[~, I_1] = min(x_stat);
x_LE_upDown = [x_LE_up x_LE_down];
[~, I_2] = min(x_LE_upDown);
deltax_shift = x_stat(I_1) - x_LE_upDown(I_2);
x_stat = x_stat - deltax_shift;


%% Shift dei profili e generazione urto pala superiore
% Shock pala inferiore e pala superiore
x_shock = [flip(x_shock_up) x_shock_down];

%% Distance between leading-edge and bow shock
deltax = min(x_stat) - min(x_shock);
x0 = [x0_down x_0_up];
SB = [x_SB_down x_SB_up; y_SB_down y_SB_up];
S = [x_S_down x_S_up; y_S_down y_S_up];
phi_up = abs(atan(gradient(y_shock_up,x_shock_up)));
phi_down = abs(atan(gradient(y_shock_down,x_shock_down)));

end


