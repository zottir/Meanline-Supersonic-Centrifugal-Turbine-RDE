function [mu_gasMixture] = gasMixtureViscosity(T,dati_mixture, dati_GasViscosityDatabase, ViscosityMixtureFormula)
% Function to calculate mixture thermal conductivity with Mason-Saxena
% model
% INPUT
% T: Temperature at which evaluate the properties [K]
% dati_mixture: table of the type:
% | Molecule   Molecular_Mass   Mass_Fraction   Molar_Fraction|
% |-----------------------------------------------------------|
% |-----------------------------------------------------------|
% dati_GasViscosityDatabase: table of the type: 
% | Fluid   mu0   T0mu   Smu   k0   T0k   Sk |
% |------------------------------------------|
% |------------------------------------------|
% ViscosityMixtureFormula: 1 Graham, 2 Wilke-Golubov, 3 Herning and Zipperer, 4 Brokaw, 5 Davidson
%
% OUTPUT
% mu_gasMixture: Mixture viscosity


%% Import mixture data
B=table2array(dati_mixture(:,2:5));
molecole=convertCharsToStrings(table2array(dati_mixture(1:end,1)))';
Mm=B(:,1)';      
X=B(:,4)';       

%% Molecule viscosity with Sutherland
mu_vec=[];
clear i
for i=1:max(size(molecole))
    gasName_i=molecole(i);
    [mu_i, k_i] = Sutherland_database(T,gasName_i, dati_GasViscosityDatabase);
    mu_vec=[mu_vec mu_i];
end

%% Graham Model
if ViscosityMixtureFormula==1
    mu_gasMixture=sum(X.*mu_vec);
end

%% Wilke-Golubov Model

if ViscosityMixtureFormula==2
    clear i
    phi=[];
    for i=1:max(size(molecole))
        for j=1:max(size(molecole))
            phi(i,j)=((1 + (((mu_vec(i)/mu_vec(j))^(1/2))*((Mm(j)/Mm(i))^(1/4))))^2)/((4/sqrt(2))*((1+(Mm(i)/Mm(j)))^(1/2)));
        end
    end 
    mu_gasMixture=sum((X.*mu_vec)./(X + (sum((((X'*ones(1, max(size(molecole))))').*(phi-diag(diag(phi)))), 2)')));
end

%% Herning-Zipperer Model
if ViscosityMixtureFormula==3
    mu_gasMixture=(sum(mu_vec.*X.*sqrt(Mm)))/(sum(X.*sqrt(Mm)));
end

%% Brokaw Model
if ViscosityMixtureFormula==4
    A_Br=[];
    m_Br=[];
    clear i j
    for i=1:max(size(molecole))
        for j=1:max(size(molecole))
            m_Br(i,j)=((4*Mm(i)*Mm(j))/((Mm(i)+Mm(j))^2))^(1/4);
            A_term=((Mm(i)/Mm(j)) - ((Mm(i)/Mm(j))^(0.45)))/((2*(1 + (Mm(i)/Mm(j)))) + (((1 + ((Mm(i)/Mm(j))^(0.45)))*m_Br(i,j))/(1+m_Br(i,j))));
            A_Br(i,j)=m_Br(i,j)*((Mm(j)/Mm(i))^(1/2))*(1 + A_term);
        end
    end
    A_Br=A_Br-diag(diag(A_Br)); % Delete diagonal elements to have i != j
    mu_gasMixture=sum((X.*sqrt(mu_vec))./((X./sqrt(mu_vec)) + (sum(((((X./sqrt(mu_vec))'*ones(1,7))').*(A_Br)), 2)')));    
end

%% Davidson Model
if ViscosityMixtureFormula==5
    fluidity_D=0;
    E_D=[];
    y_D_tot=sum(X.*sqrt(Mm));
    A_D=0.375;
    clear i j
    for i=1:max(size(molecole))
        for j=1:max(size(molecole))
            y_i_D=(X(i)*sqrt(Mm(i)))/y_D_tot;
            y_j_D=(X(j)*sqrt(Mm(j)))/y_D_tot;
            m_i_D=((X(i))^2)*Mm(i);
            m_j_D=((X(j))^2)*Mm(j);
            E_D(i,j)=(2*sqrt(m_i_D)*sqrt(m_j_D))/(m_i_D+m_j_D);
            %E_ji_D=(2*m_i_D)/(m_i_D+m_j_D);
            %E_ij_D=(2*m_j_D)/(m_i_D+m_j_D);
            %E_D_ij_alt=(y_i_D*E_ij_D + y_j_D*E_ji_D)/(y_i_D+y_j_D);
            fluidity_D=fluidity_D + ((y_i_D*y_j_D*(E_D(i,j)^A_D))/(sqrt(mu_vec(i))*sqrt(mu_vec(j))));
        end
    end
    mu_gasMixture=1/fluidity_D;   
end

end

