%% Test data base
clear all;close all;
%% This global is needed by the Cp, Cv etc functions
addpath('matlab/General/Nasa'); % Add directory of Nasa routines to Matlab-path
global Runiv Pref Tref
Runiv = 8.314;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;
kmol=1e3;
dm=0.1;
bara=1e5;
kPa = 1000;
kN=1000;
kg=1;
s=1;
cm3 = 10^(-6);
%% 
DBdir = 'General\Nasa';
DBname = 'NasaThermalDatabase';
load(fullfile(DBdir,DBname));
%%
iElements = myfind({Sp.Name},{'C2H5OH','Gasoline','H2O', 'CO2', 'N2', 'O2'});
Elements = Sp(iElements);
NElements = length(Elements);
Mi = [Elements.Mass];

Tamb = 293;                                 % [K]       Ambient temperature
Pamb = 100*kPa;                             % [Pa]      Ambient pressure

%Using the online available VOLUME ratios of air components
Vair = [0 0 0 0.000407 0.78084 0.20946];    % [V/Vtot]  Volume compisition of air
%Multiplying this by the density leaves us with mass ratios (volume * mass/volume = mass)
%However, these are not percentages that add up to 1.
%By the way, I could not find a density table in the nasa tables, so I used
%the engineering toolbox
Densities = [789 737 0.804 1.842 1.165 1.331];
YairNonunity = Vair .* Densities;   % Mass composition of air, total does NOT add up to unity
%Dividing this by the total of these ratios leaves us with a matrix in
%which every value adds up to 1. This, then, is a matrix with mass
%percentages.
Yair = YairNonunity ./ sum(YairNonunity, 'all'); % [m/mtot]   Mass composition of air (DOES add up to unity)
%Dividing this by the molar masses to obtain the mole ratios (mass / mass/moles = moles)
NairNonunity = Yair ./ Mi;
%Doing the same thing as before to make sure everything is in percentages
Nair = NairNonunity ./ sum(NairNonunity, 'all'); % [N/Ntot]  Mole composition of air

TR = [200:1:5000];

for i=1:length(Elements)
    si(:,i) = SNasa(TR, Elements(i));
    hi(:,i) = HNasa(TR, Elements(i));
    ui(:,i) = UNasa(TR, Elements(i));
end



%%
%Insert volume ratios (i.e. 5% ethanol ==> 0.05)
%'Evol' = Ethanol volume (ratio)
EvolE5          = 0.05;
EvolE10         = 0.1;
EvolE15         = 0.15;

XfuelE5         = [EvolE5 1-EvolE5 0 0 0 0]; % [V/Vtot]     Volume ratio
XfuelE10        = [EvolE10 1-EvolE10 0 0 0 0];
XfuelE15        = [EvolE15 1-EvolE15 0 0 0 0];

%See the computations at the air ratios to find out how this works.
YfuelE5Nonunity  = XfuelE5 .* Densities;
YfuelE10Nonunity = XfuelE10 .* Densities;
YfuelE15Nonunity = XfuelE15 .* Densities;

Yfuel5           = YfuelE5Nonunity ./ sum(YfuelE5Nonunity, 'all');
%   \/ I sort of hate Matlab for this reason
Yfuel10          = YfuelE10Nonunity ./ sum(YfuelE10Nonunity, 'all');
Yfuel15          = YfuelE15Nonunity ./ sum(YfuelE15Nonunity, 'all');

NfuelE5Nonunity = Yfuel5 ./ Mi;
NfuelE10Nonunity = Yfuel10 ./ Mi;
NfuelE15Nonunity = Yfuel15 ./ Mi;

NfuelE5 = NfuelE5Nonunity ./ sum(NfuelE5Nonunity, 'all');
NfuelE10 = NfuelE10Nonunity ./ sum(NfuelE10Nonunity, 'all');
NfuelE15 = NfuelE15Nonunity ./ sum(NfuelE15Nonunity, 'all');

%%
% For now using air rates but could also be simply mass in one
% thermodynamic cycle. Review later. I am using this because this is
% simalar to the Jet engine assignmed in thermo. (that is an open system.
% And this sorta is an open system but also a closed system as the air
% leaves eventually. When using the mass of fuel and air per cycle just
% have to do it times number of cycles per unit time(half the number of revolutions on
% the crank) 

%CxHyOz + (x + y / 4 - z / 2) O2 + other stuff in the air relative to the oxygen --> xCO2 + y / 2 H2O + other stuff in the air
%ALL oxygen in the air is burnt up (stoichiometric)

%ASSSUMPTION: 100% octane
x                   = 8 * NfuelE5(2) + 2 * NfuelE5(1);          %[atoms]    number of CARBON atoms in a fuel molecule
y                   = 18 * NfuelE5(2) + 6 * NfuelE5(1);         %[atoms]    number of HYDROGEN atoms in a fuel molecule
z                   = 0 * NfuelE5(2) + 1 * NfuelE5(1);          %[atoms]    number of OXYGEN atoms in a fuel molecule
molarmassFuel       = 114.232 * NfuelE5(2) + 30.07 * NfuelE5(1);%[g/mol]    G/MOL, NOT KG/MOL!!! how heavy one mole of fuel is
%airDensity          = 101325;                                  %[kg/m^3]   density of the air at the specified temperature and elevation without fuel

fuelConsumption     = 2.4333 * 10^(-4);                         %[kg/s]     how many kilogrammes of fuel are needed to let the engine operate for one second
molarRatioOxygenFuel= x + y / 4 - z / 2;                        %[atoms]    how many oxygen atoms are needed to completely combust one fuel molecule, may be a fraction.
ratioNitrogenOxygen = Yair(5) / Yair(6);                        %[ratio]    mass ratio of nitrogen to oxygen in the air 
ratioCO2Oxygen      = Yair(4) / Yair(6);                        %[ratio]    mass ratio of nitrogen to oxygen in the air 
%I neglect the other noble gasses present in the air
AFstoi              = ( molarRatioOxygenFuel * 31.9988 + ratioNitrogenOxygen * molarRatioOxygenFuel * 28.0134 + ratioCO2Oxygen * molarRatioOxygenFuel * 44.01) / (molarmassFuel * 1); 

mfurate             = fuelConsumption;                          %[kg/s]    Fuel rate required  - determine
AF                  = AFstoi;                                   %[ratio]       Air to fuel ratio - assumed to be stoichiometric
%% 1-2 Air in take - isobaric
P1 = Pamb;
T1 = Tamb; 
P2 = P1;
V2 = 196*cm3;
r = 2.5;
V1 = V2/r;

% r = V2/V1

% Using Ideal Gas Law it is possible to determine temp P1*V1/T1 = P2*V2/T2
% d = 

T2 = T1*V2/V1;  % T2 = T1 * r



%% 2-3 Adiabatic compression
mairrate            = mfurate*AF;                           % [kg/s]    Air rate coming into the carb
mrate_in            = mairrate + mfurate;                   % {kg/s]    Rate of intake air

mO2rate_in          = mairrate*Yair(6);                     % [kg/s]
mN2rate_in          = mairrate*Yair(5);                     % [kg/s]
mC2H5OHrate_in      = mfurate*Yfuel10(1);
mGASOLINErate_in    = mfurate*Yfuel10(2);
Ypre_comb           = [mC2H5OHrate_in, mGASOLINErate_in, 0, 0, mN2rate_in, mO2rate_in]./mrate_in; %mass fraction

% moles
MoleO2rate_in       = mO2rate_in/Mi(6);
MoleN2rate_in       = mN2rate_in/Mi(5);
MoleC2H5OHrate_in   = mC2H5OHrate_in/Mi(1);
MoleGASOLINErate_in = mGASOLINErate_in/Mi(2);
Molerate_in         = MoleO2rate_in + MoleN2rate_in + MoleC2H5OHrate_in + MoleGASOLINErate_in;
Xpre_comb           =[MoleGASOLINErate_in, MoleC2H5OHrate_in, 0, 0, MoleN2rate_in, MoleO2rate_in]./Molerate_in;
Mpre_comb           = Xpre_comb*Mi';
Rpre_comb           = Runiv/Mpre_comb';

V3 = V1;
Delta_V = V2 - V3;
steps =  1000;
dV = Delta_V/steps;
V23(1) = V2;
T23(1) = T1; %T2
P23(1) = P2;

for i=2:steps+1
    for ii = 1:NElements
       Cp(ii) = CpNasa(T23(i-1), Elements(ii));
       Cv(ii) = CvNasa(T23(i-1), Elements(ii));
    end
    Cp23(i) = Cp*Ypre_comb';
%     Cv23 = Cp*Xpre_comb';
    Cv23(i) = Cv*Ypre_comb';
    R(i) = Cp23(i) - Cv23(i);
    k(i) = Cp23(i)/Cv23(i);
    V23(i) = V23(i-1) - dV;
    T23(i) = T23(i-1)*(V23(i-1)/V23(i))^(k(i)-1);
    P23(i) = P23(i-1)*(V23(i-1)/V23(i))^(k(i));
%     T23(i) =  T23(i-1)*exp(Rpre_comb/Cv23(i))*(V23(i-1)/V23(i));
end

check = V23(steps+1) - V3;
T3 = T23(steps+1);
P3 = P23(steps+1)/kPa;
%%

% ds = cv ln(T2/T1) + R ln(v2/v1)

% cv(T) ln(T2/T1) =  R ln(V1/V2) --> T2 = T1*e^(R/Cv)*(V1/V2)

%s2 - s1  = 0

%s2 - sref = SNasa(T2) -Rln(P2/Pref)
%s1 - sref = SNasa(T1) -Rln(p1/Pref)


%% 3-4 Combustion - isochoric

%==========================================================================

% Q = W + dU
% For isochoric no work is being done
% Q = dU = U4-U3

%==========================================================================

for i=1:NElements
    u3(i) = UNasa(T3, Elements(i));
    h_comb(i) = HNasa(Tref, Elements(i));
end
U3 = u3*Ypre_comb'


%mass of gasoline is reported to be lighter in NASA tables then google's
%answer; Google (and my hand calc report 114 g/mol, NASA report 106.4 g/mol
%I am using the NASA one for the sake of ease

%==========================================================================

% Stoichemetric relations

% C2H5OH + 3*O2 ? 2*CO2 + 3*H2O 

% C8H18 + 12.5*O2 ? 8*CO2 +9*H2O

%==========================================================================
mH2Orate_out = 3*Mi(3)/Mi(1)*mC2H5OHrate_in + 9*Mi(3)/Mi(2)*mGASOLINErate_in;
mO2rate_out = mO2rate_in - 3*Mi(6)/Mi(1)*mC2H5OHrate_in - 12.5*Mi(6)/Mi(2)*mGASOLINErate_in;
mCO2rate_out = 2*Mi(4)/Mi(1)*mC2H5OHrate_in + 8*Mi(4)/Mi(2)*mGASOLINErate_in;
mN2rate_out = mN2rate_in;
moutrate = [0 0 mH2Orate_out mCO2rate_out, mN2rate_out, mO2rate_out];
check = sum(moutrate) - mrate_in;
Yafter_comb = moutrate/sum(moutrate);

hpre_comb = h_comb*Ypre_comb';
hafter_comb = h_comb*Yafter_comb';
Q_comb = hafter_comb - hpre_comb

U4 = U3 + Q_comb;
Uafter_comb = ui*Yafter_comb';
T4 = interp1(Uafter_comb,TR,U4)

T34(1) = T3;
P34(1) = P3*kPa;
V34(1) = V3;
for i = 2:steps+1
    for ii =1:NElements
        Cv(ii) =CvNasa(T34(i-1), Elements(ii));
    end
    Cv34(i) = Cv*Yafter_comb';
    T34(i) = -Q_comb/steps/(Cv34(i)) +T34(i-1);
    V34(i) = V3;
    P34(i) = P34(i-1)*T34(i)/T34(i-1);
end
T4= T34(steps+1)
P4 = P34(steps+1)
% p3/T3 = P4/T4

%%
V4 = V3;
V5 = V2;
Delta_V = V5 - V4;
steps =  1000;
dV = Delta_V/steps;
V45(1) = V4;
T45(1) = T4; %T2
P45(1) = P4;

for i=2:steps+1;
    for ii = 1:NElements
       Cp(ii) = CpNasa(T45(i-1), Elements(ii));
       Cv(ii) = CvNasa(T45(i-1), Elements(ii));
    end
    Cp45(i) = Cp*Yafter_comb';
%     Cv23 = Cp*Xpre_comb';
    Cv45(i) = Cv*Yafter_comb';
    R(i) = Cp45(i) - Cv45(i);
    k(i) = Cp45(i)/Cv45(i);
    V45(i) = V45(i-1) + dV;
    T45(i) = T45(i-1)*(V45(i-1)/V45(i))^(k(i)-1);
    P45(i) = P45(i-1)*(V45(i-1)/V45(i))^(k(i));
%     T23(i) =  T23(i-1)*exp(Rpre_comb/Cv23(i))*(V23(i-1)/V23(i));
end

check = V45(steps+1) - V4;
T5 = T45(steps+1);
P5 = P45(steps+1);

%%
P56(1) = P5;
Delta_P = (P5- P2);
dP = Delta_P/steps;
V56(1) = V5;
for i=2:steps
  V56(i) = V5;
  P56(i) = P56(i-1) - dP;
end

%%
hold on; plot(V23, P23); plot(V34, P34); plot(V45, P45); plot(V56, P56);
xlabel("V [m^3]");
ylabel("p [N/m^2]")
title("Theoretical p-V diagram Otto cycle assumptions")