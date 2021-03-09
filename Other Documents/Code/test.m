clear all;clc
%% add general to matlab path 
%%
run('Real_script.m');
DataDir='Data\No_fuel';ColumnOrder={'time','Sensor','Encoder'};
% DataDir='../Data/Gasoline';cOrder={'time'','Encoder','Sensor};
col = lines(3);
%% Loading all measurments in DataDir
Files=dir(fullfile(DataDir,'*.txt'));nFiles=length(Files);                  % dir gives a directory listing, only *.txt files in this case
for i=1:nFiles
    fname       = Files(i).name;                                            % Take a name from the list
    curfilename = fullfile(DataDir,fname);                                  % Create the full name
    Data        = ImportData4GB10(curfilename,ColumnOrder);                             % Read the data. Type help ImportData4GB10
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
    for ii = 2:length(V)
        Ca(ii) = Ca(ii-1) +dCa;
        Volume(ii) = Vcyl(Ca(ii));
    end
     [maxPressure, maxPressureID] = findpeaks(V, 'MinPeakHeight', 3, 'MinPeakDistance', 3000);
end

%%
% figure()
% hold on
% plot(Volume, adjustedPressure(:,1))
% plot(Volume, adjustedPressure(:,2))
% plot(Volume, adjustedPressure(:,3))
% plot(Volume, adjustedPressure(:,4))
% plot(Volume, adjustedPressure(:,5))
% plot(Volume, adjustedPressure(:,6))
% plot(Volume, adjustedPressure(:,7))
% %%
figure()
hold on
time=Data.t;
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
plot(adjustedPressure(:,1))
xlabel('Time [s]')
ylabel('Pressure [bar]')
legend('1','2','3','4','5','6','7','8','9', '10','11')

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
t=D(:,iT);p=D(:,iEncoder);V=D(:,iSensor);
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


function [V, V_Theta] = volumeCycle(Ca)
R       = 8.5;
Vt      = 196 *10^(-6);                              % [m^2]
Vc      = Vt/R;                                 % [m^2]
Vd      = Vt-Vc;
phi     =       0            %(360-153.4)/2/pi;

V=-Vd/2*cos(Ca*(2*pi/360))+Vc+Vd/2;






r = 0.030;
l = 0.085;
V_c = Vc
B = 0.06;
x = r*cos(Ca/360*2*pi) + sqrt(l^2 - r^2*(sin(Ca/360*2*pi))^2);

d_Theta = l + r - x;
V_Theta = pi*(B/2)^2*d_Theta +V_c;

end


 