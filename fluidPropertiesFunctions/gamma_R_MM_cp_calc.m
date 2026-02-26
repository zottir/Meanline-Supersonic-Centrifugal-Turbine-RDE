function [gamma, R_mixture, Mm_mixture, Cp_mixture] = gamma_R_MM_cp_calc(T,dati,MolarFractions,MolecularMasses,const_low,const_up)
% Function to calculate gamma, R, cp of a mixture at a given T
% INPUT
% T: Temperature at which evaluate the properties [K]
% MolarFractions: vector containing molar fractions of the mixture molecules 
% MolecularMasses: vector containing molar masses of the mixture molecules [kg/kmol]
% const_low,const_up: NASA polynomial coefficients. MAtrix of the type:
% |           a0 a1 a2 a3 ...|
% |molec1    -----------     |
% |molec2    -----------     |
% |molec3    -----------     |
%
% OUTPUT
% gamma: Mxture Specific Heat Capacity Ratio []
% R_mixture: Mixture Specific Gas Constant [J/kgK]
% Mm_mixture: MIxture Molecular Mass [kg/kmol]
% Cp_mixture: Mixture Specific Heat Capacity [J/kgK]

%% Mixture Properties
R_univ=8314.462618;

% Mixture Molar Mass
Mm_mixture=(sum(MolarFractions.*MolecularMasses))/sum(MolarFractions); 
R_mixture=R_univ/Mm_mixture;

% Mass Fractions
y=(MolarFractions.*MolecularMasses)/(sum(MolarFractions.*MolecularMasses));         
y_mat=y*ones(1,7);

% Conversion of coefficients to have cp in J/kgK and not J/molK
const_low_miscela=sum(const_low.*y_mat, 1);
const_up_miscela=sum(const_up.*y_mat, 1);

%% NASA Polynomilas

if T<=1000
    Cp_mixture = R_mixture*(sum(const_low_miscela(1:5).*[T^(0) T^(1) T^(2) T^(3) T^(4)]));
    H_miscela=(R_mixture*T)*(sum(const_low_miscela(1:6).*[T^(0) (T^(1))/2 (T^(2))/3 (T^(3))/4 (T^(4))/5 (T^(-1))/1]));
    S_miscela=(R_mixture)*(sum([const_low_miscela(1:5) const_low_miscela(7)].*[log(T) T (T^2)/2 (T^3)/3 (T^4)/4 (T^0)]));
else
    Cp_mixture = R_mixture*(sum(const_up_miscela(1:5).*[T^(0) T^(1) T^(2) T^(3) T^(4)]));
    H_miscela=(R_mixture*T)*(sum(const_up_miscela(1:6).*[T^(0) (T^(1))/2 (T^(2))/3 (T^(3))/4 (T^(4))/5 (T^(-1))/1]));
    S_miscela=(R_mixture)*(sum([const_up_miscela(1:5) const_up_miscela(7)].*[log(T) T (T^2)/2 (T^3)/3 (T^4)/4 (T^0)]));
end
    
gamma=(Cp_mixture)./(Cp_mixture-R_mixture);

end

