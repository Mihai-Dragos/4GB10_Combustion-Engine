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
%% PV diagrams - measured data

fname= ["E5_Full_load_1.txt","E5_Full_load_2.txt","E5_Full_load_3.txt","E5_Full_load_4.txt","E5_Full_load_5.txt","E5_Half_load_1.txt","E5_Half_load_2.txt","E5_Half_load_3.txt","E5_Half_load_4.txt","E5_Half_load_5.txt","E5_N0_load_1.txt","E5_N0_load_2.txt","E5_N0_load_3.txt","E5_N0_load_4.txt","E5_N0_load_5.txt","E15_Full_loaf_1.txt","E15_Full_loaf_2.txt","E15_Full_loaf_3.txt","E15_Full_loaf_4.txt","E15_Full_loaf_5.txt"]

for d = [1]
    clear V
    clear Pressure 
    clear Volume
    clear adjustedPressure_1
    clear adjustedPressure_2
    clear adjustedPressure_3

DataDir         = 'Data\E5';
ColumnOrder     = {'time','Sensor','Encoder'};
col = lines(3);

%% Loading all measurments in DataDir
%Uncomment for implementing all files in a folder
%Files=dir(fullfile(DataDir,'*.txt'));nFiles=length(Files);                  % dir gives a directory listing, only *.txt files in this case
%for i=1:
%fname       = Files(i).name;  
%Take a name from the list
curfilename = fullfile(DataDir,fname(d));                                  % Create the full name
Data        = ImportData4GB10(curfilename,ColumnOrder);              % Read the data. Type help ImportData4GB10
% store it for each case. Yes a struct that contains a struct and other
% things. Why not?
Case.Data     = Data;
Case.filename = fname(d);
Case.DataDir  = DataDir;
% preamble, put data in easy to use arrays
t      = Data.t;
p      = Data.pulse;
V      = ( (Data.Volt / 5 - 0.115) / 0.0154 );
RevEnd = Data.RevEnds;
NRevs  = Data.NRevs;

%Seperatating the different cycles in recorded data

peakDistance = 55; %[* 0.00001 or amount of data entries] excess difference between two peak pressures.
%This number differs for each file
%Because the sensor does not measure *exactly* a set number of cycles,
%There is an excess amount of entries that do not add up to a full cycle
%If the computation for N is used as stated below, this will mean every
%cycle is off by some amount of data entries. This number is that amount.
%Unclear? ask Vito.


N = (floor(length(V)/NRevs)) - peakDistance;
for ii = 1:floor(NRevs/2)
   LB = (ii-1)*2*N+1;
   UB = (2*N*ii);   
   Pressure(:,ii) = V([LB:UB]);
end
%Still needs testing made a small improvement to acocunt for all cases   
check = length(V) - N*NRevs;

%Checking where in the cycle the 
[~, maxPressureID] = max(Pressure(:,1)); 
if maxPressureID < RevEnd(1)
     %requires left shift
     signPhi = -1;
else
     %requires right shift
     signPhi = 1;
end

%Determining the corresponding Volume to pressure
dCa = 360/N;
Ca(1) = 0;
Volume(1) = Vcyl(Ca(1), signPhi);

for ii =2:length(Pressure)
    Ca(ii) = Ca(ii-1) + dCa;
    Volume(ii) = Vcyl(Ca(ii), signPhi);
   
end
Volume = Volume';

%Accounting the pressure for the drift
[~ , minVolumeID]=findpeaks(-Volume);

 for ii= 1: size(Pressure, 2)
  
      find_variable(ii) = mean(Pressure(find(round(Ca, 0) == 563), ii));
     [~, Diff_pressure(ii)] = max(diff(Pressure(:,ii)));
     minPressure_1(ii) = min(Pressure(minVolumeID, ii));
     minPressure_2(ii) = min(Pressure(Diff_pressure(ii), ii));
     adjustedPressure_1(:,ii) = Pressure(:,ii) - minPressure_1(ii) + 1.05;
     adjustedPressure_2(:,ii) = Pressure(:,ii) - minPressure_2(ii) + 1.05;
     adjustedPressure_3(:,ii) = Pressure(:,ii) - find_variable(ii) + 1.05;
end


%for graph titles
%cycles = size(minPressure);
%cycles = cycles(2);

%%
figure(1);
hold on;
plot(Volume, adjustedPressure_1(:,1));
plot(Volume, adjustedPressure_2(:,1));
plot(Volume, adjustedPressure_3(:,1));
legend("Method 1 (original)", "Method 2", "Method 3");
grid on;
grid minor;

xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')

end

figure(2)
hold on
time=Data.t;
plot(time([1:2*N]), adjustedPressure_1(:,1))
plot(time([1:2*N]), adjustedPressure_2(:,1))
plot(time([1:2*N]), adjustedPressure_3(:,1))

xlabel('Time [s]')
ylabel('Pressure [bar]')
legend("Method 1 (original)", "Method 2", "Method 3");

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
%V(1)=Vcyl(Ca(1), Vc, Vd); 
Volume(1) = Vcyl(Ca(1), signPhi);
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
%Cv =750;
%Rg = 300;
NSteps=NCa/dCa;

for i=2:NSteps

Ca(i)=Ca(i-1)+dCa;
%V(i)=Vcyl(Ca(i), Vc, Vd); % New volume for current crank-angle
V(i)=Vcyl(Ca(i), signPhi);
m(i)=m(i-1); % Mass is constant, valves are closed
dV=V(i)-V(i-1); % Volume change
[~,dQcom(i)] = HeatReleased(Ca(i), AF, mfurate, Yfuel, Yair, Mi, Runiv, Elements, Tref);
for ii = 1:NElements
       Cp(ii) = CpNasa(T(i-1), Elements(ii));
       Cv(ii) = CvNasa(T(i-1), Elements(ii));

end

Cv_mix(i)= Cv*Y_AF';
Cp_mix(i)= Cp*Y_AF';

Rg(i) = Cp_mix(i) - Cv_mix(i);

dT=(-dQcom(i)*dCa-p(i-1)*dV)/Cv_mix(i)/m(i-1); % 1st Law dU=dQ-pdV (closed system)


% adiabatic closed system with constant
% gas composition and constant Cv
T(i)=T(i-1)+dT;
p(i)=m(i)*Rg(i)*T(i)/V(i); % Gaslaw
end;

%% efficiency using Cp and Cv, that were calculated with non-constant
%%temperature
for i=2:NSteps
gamma(i) = Cp_mix(i)/Cv_mix(i);  %heat capacity ratio
Eff_otto(i) = 1-(1/r)^(gamma(i)-1); %otto efficiency
end;

Efficiency_average = mean(Eff_otto(i))

%% efficiency
T_cycle = mean(T); %Mean temperature during an cycle
for i = 1:NElements %Using NASA for specific heat values
    Cpi(:,i) = CpNasa(T_cycle,Elements(i));
    Cvi(:,i) = CvNasa(T_cycle,Elements(i));
end
Cp = Y_AF*Cpi';  %specific heat at constant pressure for fuel
Cv = Y_AF*Cvi'; %specific heat at constant volume for fuel
gamma = Cp/Cv  %heat capacity ratio
eff_otto = 1-(1/r)^(gamma-1) %otto efficiency
%%
%eff = trapz(dV,p)/(q_lhv*Mfuel); %Thermal efficiency

% function V = Vcyl(Ca, Vc, Vd)
% % V         - Volume at give crank angle            - [m^3]
% % Ca        - Crank angle                           - [degree]
% % Vc        - Clearance volume                      - [m^3]
% % Vd        - Displaced volume                      - [m^3]
% phi = 0;
% V=-Vd/2*cos(Ca*(2*pi/360))+Vc+Vd/2;
% 
% end

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



%%
function V = Vcyl(Ca, signPhi)
% V         - Volume at give crank angle            - [m^3]
% Ca        - Crank angle                           - [degree]
% Vc        - Clearance volume                      - [m^3]
% Vd        - Displaced volume                      - [m^3]

r       = 8.5;
Vt      = 196;                              % [m^2]
Vc      = Vt/r;                                 % [m^2]
Vd      = Vt-Vc;
phi     = signPhi*(360-153.4)/2/pi; %360-

V=-Vd/2*cos(Ca*(2*pi/360)- phi)+Vc+Vd/2;
end




% function [V, V_Theta] = volumeCycle(Ca)
% R       = 8.5;
% Vt      = 196 *10^(-6);                              % [m^2]
% Vc      = Vt/R;                                 % [m^2]
% Vd      = Vt-Vc;
% phi     = 0;            %(360-153.4)/2/pi;
% 
% V=-Vd/2*cos(Ca*(2*pi/360))+Vc+Vd/2;
% 
% r = 0.030;
% l = 0.085;
% V_c = Vc
% B = 0.06;
% x = r*cos(Ca/360*2*pi) + sqrt(l^2 - r^2*(sin(Ca/360*2*pi))^2);
% 
% d_Theta = l + r - x;
% V_Theta = pi*(B/2)^2*d_Theta +V_c;
% 
% end

function [Data]=ImportData4GB10(filename,cO)
% [Data]=ImportData4GB10(filename,cO)
%  WARNING: sometimes datafiles use comma as decimal separator. Than the
%           function FAILS. We have a fast matlab tool to convert these files
%           batchwise (see ConvertFiles.m)
%  Input 
%  [mandatory]  filename = full name to file to be read
%  [optional]   cO = cell array containing column order {'time','Encoder','Sensor'} but in
%           the order used in your data files. Above order is the default
%  Reads datafile and stores it in a struct.
%           Data.t          = time array
%           Data.pulse      = pulse array
%           Data.Volt       = sensor voltage array
%           Data.RevEnds    = Indices to revolution end point.
%           Data.NRevs      = Number of complete revolutions
% Note: function trims measurement to full cycles (so even number of revs).
%
delim='\t';
D=importdata(filename,delim, 0);
%     whos D
if (nargin < 2)
    iT=1;iEncoder=2;iSensor=3;
else
    ii = myfind(cO,{'time','Encoder','Sensor'});
    iT          =ii(1);
    iEncoder    = ii(2);
    iSensor     = ii(3);
end
t=D(:,iT);
V=D(:,iEncoder);
p=D(:,iSensor);
endloc = length(t);
% Find double tooth
[apks,alocs]    = findpeaks(p);
id              = find(apks > 0);
pks             = apks(id);
locs            = alocs(id);
tlocs           = t(locs);
dtlocs          = diff(tlocs);  % compute dt between maxima. Reason, double tooth has bigger dt!
[dlocs]         = find(dtlocs > max(dtlocs)*0.5); % Pick out double tooth
% Take full cycles (is two revolutions! this is a four-stroke engine)
CycleIndices    = locs(dlocs);
NRevs           = length(CycleIndices);
Vc              = V(CycleIndices);
if mod(NRevs,2)==0 % To make sure that only full cycles are taken.
    fprintf('%20s\n %3i revs (will be patched\n',filename,NRevs);
    if (Vc(1)<Vc(2))
        fprintf('IF condition:  %9.3f < %9.3f apparently!\n',Vc(1),Vc(2));
        CycleIndices=CycleIndices(1:end-1);
    else
        fprintf('ELSE condition:%9.3f > %9.3f apparently!\n',Vc(1),Vc(2));
        CycleIndices=CycleIndices(2:end);
    end
else
    fprintf('%20s\n %3i revs (will NOT be patched)\n',filename,NRevs);
    if (Vc(1)<Vc(2))
        fprintf('IF condition:  %9.3f < %9.3f apparently!\n',Vc(1),Vc(2));
        CycleIndices=CycleIndices(1:end-1);
    else
        fprintf('ELSE condition:%9.3f > %9.3f apparently!\n',Vc(1),Vc(2));
        CycleIndices=CycleIndices(2:end);
    end
end
%     FullCycles   = [CycleIndices(1)+1:CycleIndices(end)];
FullCycles   = [CycleIndices(1)+1:endloc];
% Store it in the struct
Data.t       = t(FullCycles)-t(FullCycles(1));  % Start at t=0. Shift it.
Data.pulse   = p(FullCycles);
Data.Volt    = V(FullCycles);
Data.RevEnds = CycleIndices(2:end)-FullCycles(1);
Data.NRevs   = length(Data.RevEnds);
end



