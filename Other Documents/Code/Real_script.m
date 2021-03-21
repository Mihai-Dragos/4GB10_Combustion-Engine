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
r       = 8.5;                                  % [-]       Compression ratio    
Vt      = 196*cm3;                              % [m^2]     Total Volume
Vc      = Vt/r;                                 % [m^2]     Clearance/dead Volume
Vd      = Vt-Vc;                                % [m^2]     Volume swept away

Xair = [0 0 0 0 0.79 0.21];                     % Molar compisition of air
Mair = Xair*Mi';                                % Molar mass of air
Yair = Xair.*Mi/Mair;                           % Mass compisition of air
mfurate = 2e-05;                                  % [kg/cyc]  Fuel per otto cycle 

% Volume compisition of fuel 
% Zfuel  =[0 1 0 0 0 0];                        % E0
% %Zfuel  =[0.05 0.95 0 0 0 0];                 % E5
% Zfuel  =[0.1 0.9 0 0 0 0];                    % E10
Zfuel  =[0.15 0.85 0 0 0 0];                    % E15

rho_C8H18 = 0.7;                                %[g/cm^3]   Density of Gasoline
rho_C2H5OH = 0.79;                              %[g/cm^3]   Density of Ethanol
rho = [rho_C2H5OH, rho_C8H18, 0 0 0 0];
Yfuel = Zfuel.*rho./(sum(Zfuel.*rho));          % Mass compistion of fuel
Mair = ((12.5*32*Yfuel(2)/114.2285) +(3*32*Yfuel(1))/46)/Yair(6);
                                                % Fuel Compisition 
Mfuel = 1;
AF = Mair/Mfuel;

Y_AF = (Yair*AF + Yfuel)/(AF+1);    % Mass composition of air-fuel mixture

TR = [200:1:5000];

for i=1:length(Elements);
    si(:,i) = SNasa(TR, Elements(i));
    hi(:,i) = HNasa(TR, Elements(i));
    ui(:,i) = UNasa(TR, Elements(i));
end

P0 = 0.537900000000000*10^5;
T0 = 293;
p(1)=0.537900000000000*10^5;
T(1)=300;
pad(1)=p(1);
Ca(1)=180;
V(1)=Vcyl(Ca(1), Vc, Vd); 
m(1) = mfurate*AF +mfurate;

%% Loop over crank-angle, with 'for' construction
NCa=360; % Number of crank-angles
dCa=0.5; % Stepsize
NSteps=NCa/dCa;
pr = 0;
Tr = 0;
Vr = 0;
p1 = p(1);
V1 = V(1);
for i=2:NSteps

Ca(i)=Ca(i-1)+dCa;
V(i)=Vcyl(Ca(i), Vc, Vd); % New volume for current crank-angle
m(i)=m(i-1); % Mass is constant, valves are closed
dV=V(i)-V(i-1); % Volume change
[~,dQcom(i), dXb(i)] = HeatReleased(Ca(i), AF, mfurate, Yfuel, Yair, Mi, Runiv, Elements, Tref);
if Ca(i) == 345
    pr = p(i-1);
    Tr = T(i-1);
    Vr = V(i-1);
  
end

for ii = 1:NElements
       Cp(ii) = CpNasa(T(i-1), Elements(ii));
       Cv(ii) = CvNasa(T(i-1), Elements(ii));
end

Cv_mix(i)= Cv*Y_AF';
Cp_mix(i)= Cp*Y_AF';
Rg(i) = Cp_mix(i) - Cv_mix(i);
k(i) = Cp_mix(i)/Cv_mix(i);
Qloss(i) = HeatLoss(Ca(i), T(i-1), p(i-1), pr, Tr, V(i-1), k(i), p1, V1, Vr);
Q = -dQcom(i)* dCa + Qloss(i)*(0.04/(720)*dCa);
dT(i)=(Q-p(i-1)*dV)/Cv_mix(i)/m(i-1); % 1st Law dU=dQ-pdV (closed system)


% adiabatic closed system with constant
% gas composition and constant Cv
T(i)=T(i-1)+dT(i);
p(i)=m(i)*Rg(i)*T(i)/V(i); % Gaslaw
end;
released = HeatReleased(Ca(i), AF, mfurate, Yfuel, Yair, Mi, Runiv, Elements, Tref)*mfurate 
figure(1)
hold on
plot(V*10^6, p/(10^5))
%% efficiency using Cp and Cv, that were calculated with non-constant
%%temperature

for i=2:NSteps
gamma(i) = Cp_mix(i)/Cv_mix(i);  %heat capacity ratio
Eff_otto(i) = 1-(1/r)^(gamma(i)-1); %otto efficiency
end;

Efficiency_average = mean(Eff_otto(i))

%% efficiency
% T_cycle = mean(T); %Mean temperature during an cycle
% for i = 1:NElements %Using NASA for specific heat values
%     Cpi(:,i) = CpNasa(T_cycle,Elements(i));
%     Cvi(:,i) = CvNasa(T_cycle,Elements(i));
% end
% Cp = Y_AF*Cpi';  %specific heat at constant pressure for fuel
% Cv = Y_AF*Cvi'; %specific heat at constant volume for fuel
% gamma = Cp/Cv  %heat capacity ratio
% eff_otto = 1-(1/r)^(gamma-1) %otto efficiency
%%
%eff = trapz(dV,p)/(q_lhv*Mfuel); %Thermal efficiency




function V = Vcyl(Ca, Vc, Vd)
% V         - Volume at give crank angle            - [m^3]
% Ca        - Crank angle                           - [degree]
% Vc        - Clearance volume                      - [m^3]
% Vd        - Displaced volume                      - [m^3]
phi = 0;
V=-Vd/2*cos(Ca*(2*pi/360))+Vc+Vd/2;

end

function [q_lhv,Qcomb,dXb] = HeatReleased(Ca, AF, mfurate, Yfuel, Yair, Mi, Runiv, Elements, Tref)

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

% hpre_comb = h_comb*Ypre_comb';
% hafter_comb = h_comb*Yafter_comb';
% q_lhv = hafter_comb - hpre_comb;


%TEST
E=15
d_octane =  703;                                                           %Density of octane [kg/m^3]
d_ethanol = 789 ;                                                          %Density of ethanol [kg/m^3]
d_fuel = (E/100)*d_ethanol + (1-E/100)*d_octane ;                          %Density of fuel [kg/m^3]

Qlhv_gasoline_1 = 44.4e6;                                                  %Lower heating value gasoline [J/kg]
Qlhv_ethanol_1 = 26.7e6;                                                   %Lower heating value ethanol [J/kg

Qlhv_gasoline_2 = Qlhv_gasoline_1  * d_octane;                             %Lower heating value gasoline [J/m^3]
Qlhv_ethanol_2 = Qlhv_ethanol_1 * d_ethanol;                               %Lower heating value ethanol [J/m^3]
Qlhv_fuel_1 =  (1-E/100)* Qlhv_gasoline_2 +(E/100)*Qlhv_ethanol_2;         %Lower heating value fuel [J/m^3]
q_lhv = -Qlhv_fuel_1/d_fuel; 



n =3;
a = 5;
Theta_d = 35;
Theta_s = (360-15);

if Ca >= Theta_s
    Xb = 1 - exp(-a*((Ca-Theta_s)/Theta_d)^n);
    dXb = n*a*(1-Xb)/Theta_d*((Ca-Theta_s)/Theta_d)^(n-1);
    dQcomb_dTheta = q_lhv*mfurate*dXb;
    Qcomb = dQcomb_dTheta*1;
elseif Ca < Theta_s  
    Qcomb =0;
    Xb = 0;
    dXb = 0;
end 
end

function Qloss = HeatLoss(Ca, T, p, pr, Tr, V, k, p1, V1, Vr)
%% variables
% Ca    - Crank angle [deg]
% Sp    - Average Piston Speed [m/s]
% Tr    - temp. start of combustion [K]
% T     - Instantaneous temp. [K]
% Tw    - Temperature of the wall [K]
% pr    - pressure start of combustion [Pa]
% p1    - Pressure start of adiabatic compression [Pa]
% V1    - Volume start of adiabatic compression [m^3]
% Vr    - Volume start of combustion [m^3]
% pm    - motored cylinder pressure at the same crank angle as 'p' [Pa]
% p     - Instantaneous pressure [pa] 
% S     - Lenght of the stroke [m]
% B     - Bore diameter [m]
% RpS   - Rotations per second [deg/s]
% k     - Heat capacity ratio (Cp/Cv) [-]
% C1,C2 - Constants for in otto cycle
if pr == 0
    Qloss =0;
else
    if Ca >= 345 && Ca <= 540
        C1 = 2.28;
        C2 = 3.24e-03;
    elseif Ca > 180 && Ca <345
        C1 = 2.28; 
        C2 =0;
    else
        C1 = 0;
        C2 = 0;
    end



    % if Ca == INSERT 
    %     pr = p
    %     Tr = T
    % end

    S   = 0.055;        %[m]
    RpS = 50;
    Sp  = 2*S*RpS;      %[m/s]
    B   = 0.067;        %[m]
    Tw  = 273.15+100;    %%% CHECK %%%
    %% Computations

    % omega - The flame speed [m/s]
    % hc    - Heat convection coefficient [w/m^2/K]
    % h     - Height of the cylinder exposed [m}
    % A     - Area of the cylinder exposed [m^2]
    % Qloss - Power lost to the wall [W, J/s]
    % pm    - Pressure of motored engine [Pa]
    Vt  = 196*10^(-6);
    Vc  = Vt/8.5;
    Vd  = Vt-Vc;

    pm = p1*(V1/V)^k;
    omega   = C1*Sp + C2*(Vd*Tr)/(pr*Vr)*(p - pm);
    hc      = 3.26*B^(-0.2)*(p/1000)^(0.8)*T^(-0.55)*omega^(0.8);
    h       = (V - Vc)/(Vd)*S;
    A       = pi/2*B^2 + 2*pi*B*h;
    if T < Tw
        Qloss = 0;
    else 
        Qloss   = hc*A*(T-Tw);
    end
    end
end




