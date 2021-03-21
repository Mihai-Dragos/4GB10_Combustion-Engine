%% Initialisation:
%  Please refrain from touching this code unless errors dictate this
clear all;                              % Clears all previously used variables
addpath('matlab/General/Nasa');         % Add directory of Nasa routines to Matlab-path
DBdir = 'General\Nasa';
DBname = 'NasaThermalDatabase';
load(fullfile(DBdir,DBname));           % Loading Nasa database

%% Global variables and units
%  Should not need to be changed
global Runiv Pref Tref kJ kmol dm bar kPa kN kg s mm K cm3
Runiv   = 8.314;     % gas constant
Pref    = 1.01235e5; % Reference pressure, 1 atm!
Tref    = 298.15;    % Reference Temperature
kJ      = 1e3;      % kilo joule
kmol    = 1e3;      % kilo mole
dm      = 0.1;      % decimetre
bar     = 1e5;      % atmospheric pressure in pascal
kPa     = 1000;     % kilo pascal
kN      = 1000;     % kilo newton
kg      = 1;
s       = 1;
mm      = 0.001;    
K       = 1;
cm3     = 10^(-6);  % cubic centimetre

%% Settings
% #########################################################################
% ######                FILL THESE IN TO CHANGE SETTINGS              #####
% #########################################################################

Tambient           = 293        *K;      % Ambient temperature
Pambient           = 100        *kPa;    % Ambient pressure
CompressionRatio   = 8.5;                % Compression ratio    
Vtotal             = 196        *cm3;    % Total Volume
EthanolPercentage  = 5          /100;    % 5 for E5, 0 for E0 ....
FuelPerCycle       = 1.73e-05   *kg;     % Fuel per otto cycle 

%These are MASS fractions, not VOLUME fractions, so 78% N2 and 21% O2 would
%be FALSE!!!
%https://www.google.com/url?sa=i&url=https%3A%2F%2Fwps.prenhall.com%2Fwps%2Fmedia%2Fobjects%2F4678%2F4790892%2Fch09_01.htm&psig=AOvVaw3aGWrax-e16Xde9JTambientik8V&ust=1616424082852000&source=images&cd=vfe&ved=0CAIQjRxqFwoTCOiKlIzPwe8CFQAAAAAdAAAAABAU
MassAir              = [0 0 0 0.0005 0.7552 0.2314];
VolumeAir            = [0 0 0 0.00037 0.7808 0.2095];
density_C8H18        = 0.7;                %[g/cm^3]   Density of Gasoline
density_C2H5OH       = 0.79;               %[g/cm^3]   Density of Ethanol
MolarMass_C8H18      = 114.2285;
MolarMass_C2H5OH     = 46.07;
ReactionCoefficients = [3*32 12.5*32 0 0 0 0];   % [Ethanol Gasoline] 

NCa             = 360;            % Number of crank-angles
stepCa          = 0.5;            % Stepsize of the loop
CaMaximum       = 345;            % Combustion angle

%% Chemical composition computations
%Air computations
iElements       = myfind({Sp.Name},{'C2H5OH', 'Gasoline', 'H2O', 'CO2', 'N2', 'O2'});
Elements        = Sp(iElements);
NElements       = length(Elements);            % Number of elements ocnsidered
Mi              = [Elements.Mass];             % Molar mass of each component
Vclearance      = Vtotal / CompressionRatio;   % [m^3]     Clearance/dead volume
Vdisplacement   = Vtotal - Vclearance;         % [m^3]     Displacement volume
MolarAir        = MassAir ./ Mi;               % Moles of air

%Fuel Computations
VolumeFuel      = [EthanolPercentage (1 - EthanolPercentage) 0 0 0 0];  % Volume compisition of fuel 
FuelDensities   = [density_C2H5OH density_C8H18 0 0 0 0];               % Densities of each component in the fuel
MassFuel        = VolumeFuel .* FuelDensities ./ ( sum(VolumeFuel .* FuelDensities) );          % Mass composition of fuel
MassAirPerMassFuel = ( (ReactionCoefficients(2) * MassFuel(2) / MolarMass_C8H18) +... 
                       (ReactionCoefficients(1) * MassFuel(1) / MolarMass_C2H5OH) )... 
                       / MassAir(6);   %How many kgs of air are needed to burn one kg of fuel optimally
AF = MassAirPerMassFuel;                               % Air to fuel ratio;
MassCompositionAF = (MassAir*AF + MassFuel)/(AF+1);    % Mass composition of air-fuel mixture

TR = [200:1:5000];

%% Iterative computations

%Computations for the enthalpy, entropy and internal energy in each state.
%Thus far, it remains unused.
for i=1:length(Elements)
    si(:,i) = SNasa(TR, Elements(i));
    hi(:,i) = HNasa(TR, Elements(i));
    ui(:,i) = UNasa(TR, Elements(i));
end

%Initial state
p(1) = 0.537900000000000*10^5;
T(1) = 300;
Ca(1)= 180;
V(1) = volumeCycle(Ca(1), Vtotal, CompressionRatio, 0); 
m(1) = FuelPerCycle*AF +FuelPerCycle;

%% Iterative computations

NSteps  = NCa/stepCa; %Number of steps made in total
pr      = 0; %These variables are empty but need to be declared outside of the for loop
Tr      = 0;
Vr      = 0;
for i=2:NSteps

    Ca(i)   = Ca(i-1) + stepCa;                                % Crank angle
    V(i)    = volumeCycle(Ca(i), Vtotal, CompressionRatio, 0); % New volume for current crank-angle
    m(i)    = m(i-1);                                          % Mass is constant, valves are closed
    dV      = V(i) - V(i-1);                                   % Volume change
    
    %See HeatReleased function
    [~,dQcombustion(i), dXb(i)] = HeatReleased(Ca(i), AF, FuelPerCycle, MassFuel, MassAir, Mi, Runiv, Elements, Tref);
    
    %When the crank angle is at its maximum, fill in the optimum values
    if Ca(i) == CaMaximum
        pr = p(i-1);
        Tr = T(i-1);
        Vr = V(i-1);
    end

    %Compute the heat capacities for each element in the mixture
    for ii = 1:NElements
        Cp(ii) = CpNasa(T(i-1), Elements(ii));
        Cv(ii) = CvNasa(T(i-1), Elements(ii));
    end

    Cv_mix(i)   = Cv * MassCompositionAF';
    Cp_mix(i) 	= Cp * MassCompositionAF';
    Rg(i)       = Cp_mix(i) - Cv_mix(i);        % Rg = cp - cv
    k(i)        = Cp_mix(i) / Cv_mix(i);        % k = cp / cv (sometimes called gamma)
    Qloss(i)    = HeatLoss(Ca(i), T(i-1), p(i-1), pr, Tr, V(i-1), k(i), p(1), V(1), Vr);
    
    % Q = -(dQcombustion per step) - Qloss of this step * (0.04(?) / 360)
    Q           = -dQcombustion(i) * stepCa - Qloss(i) * (0.04 / (2 * NCa) * stepCa);
    
    % 1st Law dU=dQ-pdV (closed system)
    % dT = (Q - pdV)/(cv * m)
    dT(i)       = (Q - p(i-1) * dV) / Cv_mix(i) / m(i-1); % 
    % T now = T previously + T change this cycle
    T(i)=T(i-1)+dT(i);
    
    % adiabatic closed system with constant
    % gas composition and constant Cv
    % p = m * Rg * T / V
    p(i)= m(i) * Rg(i) * T(i) / V(i); % Gaslaw
    
    %otto efficiency
    OttoEfficiency(i) = 1 - (1 / CompressionRatio) ^ (k(i) - 1); 
    
end

HeatRelease = HeatReleased(Ca(i), AF, FuelPerCycle, MassFuel, MassAir, Mi, Runiv, Elements, Tref) * FuelPerCycle;
Efficiency_average = mean(OttoEfficiency(i));

figure(1)
hold on
plot(V*10^6, p/(10^5))

%% Functions

function V = Vcyl(Ca, Vc, Vd)
% V         - Volume at give crank angle            - [m^3]
% Ca        - Crank angle                           - [degree]
% Vc        - Clearance volume                      - [m^3]
% Vd        - Displaced volume                      - [m^3]
phi = 0;
V   = -Vd / 2 * cos(Ca * (2 * pi / 360)) + Vc + Vd / 2;

end

function [q_lhv,Qcomb,dXb] = HeatReleased(Ca, AF, FuelPerCycle, MassFuel, MassAir, Mi, Runiv, Elements, Tref)

% Qcomb     - Energy released during combustion

AirPerCycle             = FuelPerCycle * AF;          % [kg/s]    Air rate coming into the carb
MassFlowPerCycle        = AirPerCycle + FuelPerCycle; % {kg/s]    Rate of intake air

%Mass pre-combustion
OxygenInflow            = AirPerCycle  * MassAir(6);     % [kg/s]
NitrogenInflow          = AirPerCycle  * MassAir(5);     % [kg/s]
EthanolInflow           = FuelPerCycle * MassFuel(1);
GasolineInflow          = FuelPerCycle * MassFuel(2);
MassFlowPreCombustion   = [EthanolInflow, GasolineInflow, 0, 0, NitrogenInflow, OxygenInflow] ./ MassFlowPerCycle; %mass fraction

% moles in
MoleO2rate_in           = OxygenInflow    / Mi(6);
MoleN2rate_in           = NitrogenInflow  / Mi(5);
MoleC2H5OHrate_in       = EthanolInflow   / Mi(1);
MoleGASOLINErate_in     = GasolineInflow  / Mi(2);
MoleInflow              = MoleO2rate_in + MoleN2rate_in + MoleC2H5OHrate_in + MoleGASOLINErate_in;

MoleFlowPreCombustion   = [MoleGASOLINErate_in, MoleC2H5OHrate_in, 0, 0, MoleN2rate_in, MoleO2rate_in] ./ MassFlowPerCycle;
MolarMassPreCombustion  = MoleFlowPreCombustion*Mi';

GasConstantPreCombustion= Runiv/MolarMassPreCombustion';

% after combustion;
% mass out

% 3 ethanol + 9 gasoline = water
WaterOutflow            =                3 * Mi(3)/Mi(1) * EthanolInflow + 9    * Mi(3)/Mi(2) * GasolineInflow;
% initial oxygen - 3 ethanol - 12,5 gasoline = oxygen
OxygenOutflow           = OxygenInflow - 3 * Mi(6)/Mi(1) * EthanolInflow - 12.5 * Mi(6)/Mi(2) * GasolineInflow;
% 2 ethanol + 8 gasoline = co2
CarbondioxideOutflow    =                2 * Mi(4)/Mi(1) * EthanolInflow + 8    * Mi(4)/Mi(2) * GasolineInflow;
% N2 in = N2 out
NitrogenOutflow         = NitrogenInflow;

MassOutflow             = [0 0 WaterOutflow CarbondioxideOutflow, NitrogenOutflow, OxygenOutflow];
check                   = sum(MassOutflow) - MassFlowPerCycle; %SHOULD BE ZERO
MolarMassPostCombustion = MassOutflow/sum(MassOutflow);

%moles out
MoleO2rate_out          = OxygenOutflow/Mi(6);
MoleN2rate_out          = NitrogenOutflow/Mi(5);
MoleCO2rate_out         = CarbondioxideOutflow/Mi(4);
MoleH20rate_out         = WaterOutflow/Mi(3);
MoleOutflow             = MoleO2rate_out + MoleN2rate_out + MoleCO2rate_out + MoleH20rate_out;

MoleFlowPostCombustion  = [0, 0, MoleH20rate_out, MoleCO2rate_out, MoleN2rate_out, MoleO2rate_out]./MoleOutflow;
MolarMassPostCombustion = MoleFlowPostCombustion*Mi';
GasConstantPostCombustion = Runiv/MolarMassPostCombustion';

for i=1:length(Elements)
    EnthalpyCombustion(i) = HNasa(Tref, Elements(i));
end

EnthalpyPreCombustion   = EnthalpyCombustion * MassFlowPreCombustion';
EnthalpyPostCombustion  = EnthalpyCombustion * MolarMassPostCombustion';
q_lhv                   = -EnthalpyPostCombustion + EnthalpyPreCombustion;

n = 3;
a = 5;
Theta_d = 35;
Theta_s = (360-10);

if Ca >= Theta_s
    Xb              = 1 - exp(-a * ( (Ca-Theta_s) / Theta_d ) ^ n );
    dXb             = n * a * (1-Xb) / Theta_d * ( (Ca-Theta_s) / Theta_d ) ^ (n-1);
    dQcomb_dTheta   = q_lhv * FuelPerCycle * dXb;
    Qcomb           = dQcomb_dTheta*1;
elseif Ca < Theta_s  
    Qcomb           = 0;
    Xb              = 0;
    dXb             = 0;
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
    Tw  = 273.15+95;    %%% CHECK %%%
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
        Qloss   = -hc*A*(T-Tw);
    end
    end
end