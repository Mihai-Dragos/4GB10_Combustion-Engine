% Initialisation
%% Test data base
clear all;
% close all;
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
mm  = 0.001;    %[m]
cm3 = 10^(-6);  %[m^2]

%% 
DBdir = 'General\Nasa';
DBname = 'NasaThermalDatabase';
load(fullfile(DBdir,DBname));
%% Constants
iElements = myfind({Sp.Name},{'C2H5OH','Gasoline','H2O', 'CO2', 'N2', 'O2'});
Elements = Sp(iElements);
NElements = length(Elements);
Mi = [Elements.Mass];

Tamb    = 293;                                  % [K]       Ambient temperature
Pamb    = 100*kPa;                              % [Pa]      Ambient pressure
r       = 8.5;
Vt      = 196*cm3;                              % [m^2]
Vc      = Vt/r;                                 % [m^2]
Vd      = Vt-Vc;

Xair = [0 0 0 0 0.79 0.21];                 % Molar compisition of air
Mair = Xair*Mi';                            % Molar mass of air
Yair = Xair.*Mi/Mair;                       % Mass compisition of air

mfurate = 0.1;

%Zfuel  =[0 1 0 0 0 0];             %E0
%Zfuel  =[0.05 0.95 0 0 0 0];       %E5
%Zfuel  =[0.1 0.9 0 0 0 0];          %E10
Zfuel  =[0.15 0.85 0 0 0 0];       %E15

rho_C8H18 = 0.7; %[g/cm^3];
rho_C2H5OH = 0.79; %[g/cm^3]
rho = [rho_C2H5OH, rho_C8H18, 0 0 0 0];
Yfuel = Zfuel.*rho./(sum(Zfuel.*rho));
Mair = ((12.5*32*Yfuel(2)/114.2285) +(3*32*Yfuel(1))/46)/Yair(6);
Mfuel = 1;
AF = Mair/Mfuel;

Y_AF = (Yair*AF + Yfuel)/(AF+1);    % Mass composition of air-fuel mixture

TR = [200:1:5000];

for i=1:length(Elements);
    si(:,i) = SNasa(TR, Elements(i));
    hi(:,i) = HNasa(TR, Elements(i));
    ui(:,i) = UNasa(TR, Elements(i));
end

P0 = Pamb;
T0 = Tamb;
p(1)=P0;
T(1)=T0;
pad(1)=p(1);
Ca(1)=180;
V(1)=Vcyl(Ca(1), Vc, Vd); 
m(1) = 1;


%% Vcyl is a function that
% computes cyl vol as fie of crank-
% angle for given B,S,l and rc
%m(1)=p(1)*V(1)/Rg/T(1);
for i = 340:400
[~,Qcomb(i)] = HeatReleased(i, AF, mfurate, Yfuel, Yair, Mi, Runiv, Elements, Tref);
end
%% Loop over crank-angle, with 'for' construction
NCa=360; % Number of crank-angles
dCa=0.5; % Stepsize
Cv =750;
Rg = 300;
NSteps=NCa/dCa;
for i=2:NSteps
Ca(i)=Ca(i-1)+dCa;
V(i)=Vcyl(Ca(i), Vc, Vd); % New volume for current crank-angle
m(i)=m(i-1); % Mass is constant, valves are closed
dV=V(i)-V(i-1); % Volume change
[~,dQcom(i)] = HeatReleased(Ca(i), AF, mfurate, Yfuel, Yair, Mi, Runiv, Elements, Tref);
dT=(-dQcom(i)*dCa-p(i-1)*dV)/Cv/m(i-1); % 1st Law dU=dQ-pdV (closed system)
% adiabatic closed system with constant
% gas composition and constant Cv
T(i)=T(i-1)+dT;
p(i)=m(i)*Rg*T(i)/V(i); % Gaslaw
end;

% efficiency
T_cycle = mean(T); %Mean temperature during an cycle
for i = 1:NElements %Using NASA for specific heat values
    Cpi(:,i) = CpNasa(T_cycle,Elements(i));
    Cvi(:,i) = CvNasa(T_cycle,Elements(i));
end
Cp = Y_AF*Cpi';  %specific heat at constant pressure for fuel
Cv = Y_AF*Cvi'; %specific heat at constant volume for fuel
gamma = Cp/Cv  %heat capacity ratio
eff_otto = 1-(1/r)^(gamma-1) %otto efficiency

%eff = trapz(dV,p)/(q_lhv*Mfuel); %Thermal efficiency


function V = Vcyl(Ca, Vc, Vd)
% V         - Volume at give crank angle            - [m^3]
% Ca        - Crank angle                           - [degree]
% Vc        - Clearance volume                      - [m^3]
% Vd        - Displaced volume                      - [m^3]
phi = 0;
V=-Vd/2*cos(Ca*(2*pi/360))+Vc+Vd/2;

end

function [q_lhv,Qcomb] = HeatReleased(Ca, AF, mfurate, Yfuel, Yair, Mi, Runiv, Elements, Tref)

% Qcomb     - Energy released during combustion

mairrate      	 = mfurate*AF;               % [kg/s]    Air rate coming into the carb
mrate_in         = mairrate + mfurate;       % {kg/s]    Rate of intake air

%Mass pre-combustion
mO2rate_in              = mairrate*Yair(6);     % [kg/s]
mN2rate_in              = mairrate*Yair(5);     % [kg/s]
mC2H5OHrate_in          = mfurate*Yfuel(1);
mGASOLINErate_in        = mfurate*Yfuel(2);
Ypre_comb               = [mC2H5OHrate_in, mGASOLINErate_in, 0, 0, mN2rate_in, mO2rate_in]./mrate_in; %mass fraction

% moles in
MoleO2rate_in = mO2rate_in/Mi(6);
MoleN2rate_in = mN2rate_in/Mi(5);
MoleC2H5OHrate_in = mC2H5OHrate_in/Mi(1);
MoleGASOLINErate_in = mGASOLINErate_in/Mi(2);
Molerate_in = MoleO2rate_in + MoleN2rate_in + MoleC2H5OHrate_in + MoleGASOLINErate_in;
Xpre_comb =[MoleGASOLINErate_in, MoleC2H5OHrate_in, 0, 0, MoleN2rate_in, MoleO2rate_in]./Molerate_in;
Mpre_comb = Xpre_comb*Mi';

Rpre_comb = Runiv/Mpre_comb';

% after combustion;
%mass out
mH2Orate_out = 3*Mi(3)/Mi(1)*mC2H5OHrate_in + 9*Mi(3)/Mi(2)*mGASOLINErate_in;
mO2rate_out = mO2rate_in - 3*Mi(6)/Mi(1)*mC2H5OHrate_in - 12.5*Mi(6)/Mi(2)*mGASOLINErate_in;
mCO2rate_out = 2*Mi(4)/Mi(1)*mC2H5OHrate_in + 8*Mi(4)/Mi(2)*mGASOLINErate_in;
mN2rate_out = mN2rate_in;
moutrate = [0 0 mH2Orate_out mCO2rate_out, mN2rate_out, mO2rate_out];
check = sum(moutrate) - mrate_in;
Yafter_comb = moutrate/sum(moutrate);

%moles out
MoleO2rate_out = mO2rate_out/Mi(6);
MoleN2rate_out = mN2rate_out/Mi(5);
MoleCO2rate_out = mCO2rate_out/Mi(4);
MoleH20rate_out = mH2Orate_out/Mi(3);
Molerate_out = MoleO2rate_out + MoleN2rate_out + MoleCO2rate_out + MoleH20rate_out;
Xaft_comb =[0, 0, MoleH20rate_out, MoleCO2rate_out, MoleN2rate_out, MoleO2rate_out]./Molerate_out;
Maft_comb = Xaft_comb*Mi';
Raft_comb = Runiv/Maft_comb';

for i=1:length(Elements)
    h_comb(i) = HNasa(Tref, Elements(i));
end

hpre_comb = h_comb*Ypre_comb';
hafter_comb = h_comb*Yafter_comb';
q_lhv = hafter_comb - hpre_comb;

n =3;
a = 5;
Theta_d = 35;
Theta_s = (360-15);



if Ca >= Theta_s
    Xb = 1 - exp(-a*((Ca-Theta_s)/Theta_d)^n);
    dQcomb_dTheta = q_lhv*mfurate*n*a*(1-Xb)/Theta_d*((Ca-Theta_s)/Theta_d)^(n-1);
    Qcomb = dQcomb_dTheta*1;
elseif Ca < Theta_s  
    Qcomb =0;
end 
end



