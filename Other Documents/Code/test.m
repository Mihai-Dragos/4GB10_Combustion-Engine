clear all;clc
%% add general to matlab path 
addpath('General');
%%
DataDir='C:\Users\20190762\Documents\Universiteit\2020-2021 (Q3) Combustion Engine (4GB10)\4GB10_matlab_and_trainingdata\matlab\Data\No_fuel';ColumnOrder={'time','Sensor','Encoder'};
% DataDir='../Data/Gasoline';cOrder={'time'','Encoder','Sensor};
col = lines(3);
%% Loading all measurments in DataDir
Files=dir(fullfile(DataDir,'*.txt'));nFiles=length(Files);                  % dir gives a directory listing, only *.txt files in this case

for i=1:nFiles
    fname       = Files(i).name;                                            % Take a name from the list
    curfilename = fullfile(DataDir,fname);                                  % Create the full name
    [Data,locs, dtlocs]        = ImportData4GB10(curfilename,ColumnOrder);                             % Read the data. Type help ImportData4GB10
    % store it for each case. Yes a struct that contains a struct and other
    % things. Why not?
    Case(i).Data     = Data;
    Case(i).filename = fname;
    Case(i).DataDir  = DataDir;
    
    % preamble, put data in easy to use arrays
    t      = Data.t;
    p      = Data.pulse;
    V      = (((Data.Volt/5)-0.115)/0.00385);
    RevEnd = Data.RevEnds;
    NRevs  = Data.NRevs;
   
    %Seperatating the different cycles in recorded data
    N = (floor(length(V)/NRevs));
    
    
    TotalRotation = NRevs*360;
    dCa = TotalRotation/length(V);
    Ca(1) = 0;
    for ii = 2:length(V);
        Ca(ii) = Ca(ii-1) +dCa;
        Volume(ii) = Vcyl(Ca(ii));
    end
     [maxPressure, maxPressureID] = findpeaks(V, 'MinPeakHeight', 3, 'MinPeakDistance', 3000);
     
end

%%
figure()
hold on
plot(Volume, adjustedPressure(:,1))
plot(Volume, adjustedPressure(:,2))
plot(Volume, adjustedPressure(:,3))
plot(Volume, adjustedPressure(:,4))
plot(Volume, adjustedPressure(:,5))
plot(Volume, adjustedPressure(:,6))
plot(Volume, adjustedPressure(:,7))
%%
figure()
hold on
time=Data.t
plot(time([1:2*N]), Pressure(:,1))
plot(time([1:2*N]), Pressure(:,2))
plot(time([1:2*N]), Pressure(:,3))
plot(time([1:2*N]), Pressure(:,4))
plot(time([1:2*N]), Pressure(:,5))
plot(time([1:2*N]), Pressure(:,6))
plot(time([1:2*N]), Pressure(:,7))
plot(time([1:2*N]), Pressure(:,8))
plot(time([1:2*N]), Pressure(:,9))
plot(time([1:2*N]), Pressure(:,10))
plot(time([1:2*N]), Pressure(:,11))
% plot(adjustedPressure(:,))
% xlabel('Time [s]')
% ylabel('Pressure [bar]')
% legend('1','2','3','4','5','6','7','8','9', '10','11')

%%
function V = Vcyl(Ca)
% V         - Volume at give crank angle            - [m^3]
% Ca        - Crank angle                           - [degree]
% Vc        - Clearance volume                      - [m^3]
% Vd        - Displaced volume                      - [m^3]

r       = 8.5;
Vt      = 196;                              % [m^2]
Vc      = Vt/r;                                 % [m^2]
Vd      = Vt-Vc;
phi     = -1*(360-153.4)/2/pi;

V=-Vd/2*cos(Ca*(2*pi/360)- phi)+Vc+Vd/2;
end
 