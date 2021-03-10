clear all;
clc

%% add general to matlab path 
%%
%run('Real_script.m');
fname= ["E5_Full_load_1.txt","E5_Full_load_2.txt","E5_Full_load_3.txt","E5_Full_load_4.txt","E5_Full_load_5.txt","E5_Half_load_1.txt","E5_Half_load_2.txt","E5_Half_load_3.txt","E5_Half_load_4.txt","E5_Half_load_5.txt","E15_Full_loaf_1.txt","E15_Full_loaf_2.txt","E15_Full_loaf_3.txt","E15_Full_loaf_4.txt","E15_Full_loaf_5.txt","E15_Full_loaf_6.txt","E15_half_load_1.txt","E15_half_load_2.txt","E15_half_load_3.txt","E15_half_load_4.txt","E15_half_load_5.txt"];
% fname(1)="E5_Full_load_1.txt";
% fname(2)="E5_Half_load_1.txt";
%file=['Data\E5','Data\E10']
d=0;
for d=[1 2 3 4 6 7 8 9 10]
    clear V
    clear Pressure 
    clear Volume
    clear adjustedPressure
%  if d <11
 fname(d);
 DataDir         = 'Data\E5';
%  else
%   DataDir = 'Data\E15';
%   end
 

% fname           = "E5_Full_load_1.txt";
ColumnOrder     = {'time','Sensor','Encoder'};
% DataDir='../Data/Gasoline';cOrder={'time'','Encoder','Sensor};
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
V      = ( (Data.Volt / 5 - 0.115) / 0.00385 );
RevEnd = Data.RevEnds;
NRevs  = Data.NRevs;

%finding the double tooth position
if Data.NRevs >= 2
    [~, maxPressureID] = findpeaks(Data.Volt, 'MinPeakHeight', 0.5, 'MinPeakDistance', 3000);
    stepAngle = 360/(abs(Data.RevEnds(1) - Data.RevEnds(2)));
    if 0 == isempty(maxPressureID)
        for ii = 1:2
           stepDiff=abs(maxPressureID(2) - Data.RevEnds(ii));
           angle(ii) = stepAngle * stepDiff;
        end
    end
    meanAngle = mean(angle);
    stdAngle = std(angle);
end



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
if maxPressureID < RevEnd(2*(1-1)+1)
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
[~, minVolumeID]=findpeaks(-Volume);
 for ii = 1:size(Pressure, 2)
     minPressure(ii) = min(Pressure(minVolumeID, ii));
     adjustedPressure(:,ii) = Pressure(:,ii) - minPressure(ii) + 1.05;
 end

%Shifting the peak pressure


% Actual Plotting use of subplot
% figure(1)
% subplot(1,2,2);                                          % help subplot if you don't understand
% plot(t,p,'b-');
% line(t(RevEnd),p(RevEnd),'LineStyle','none','Marker','s','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','y');
% xlabel('t [s]');ylabel('pulse [V]');
% legend('pulse','double tooth');
% title(fname);
% 
% subplot(1,2,1);
% plot(t,V,'r-');
% line(t(RevEnd),V(RevEnd),'LineStyle','none','Marker','s','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','y');
% xlabel('t [s]');ylabel('p Sensor Signal [V]');
% legend('p signal','location of double tooth');
% title(fname);

%for graph titles
cycles = size(minPressure);
cycles = cycles(2);

%%
figure(1)
hold on
plot(Volume, adjustedPressure(:,1))
% plot(Volume, adjustedPressure(:,2))
% plot(Volume, adjustedPressure(:,3))
% plot(Volume, adjustedPressure(:,4))
% plot(Volume, adjustedPressure(:,5))
% plot(Volume, adjustedPressure(:,6))
% plot(Volume, adjustedPressure(:,7))
% plot(Volume, adjustedPressure(:,8))
% plot(Volume, adjustedPressure(:,9))
% plot(Volume, adjustedPressure(:,10))
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')
%title(fname + " (" + cycles + " cycles)");
legend(fname)
%%
figure(2)
hold on
time=Data.t;
plot(time([1:2*N]), adjustedPressure(:,1))
% plot(time([1:2*N]), Pressure(:,2))
% plot(time([1:2*N]), Pressure(:,3))
% plot(time([1:2*N]), Pressure(:,4))
% plot(time([1:2*N]), Pressure(:,5))
% plot(time([1:2*N]), Pressure(:,6))
% plot(time([1:2*N]), Pressure(:,7))
% plot(time([1:2*N]), Pressure(:,8))
% plot(time([1:2*N]), Pressure(:,9))
% plot(time([1:2*N]), Pressure(:,10))
% plot(time([1:2*N]), Pressure(:,11))
% plot(adjustedPressure(:,))
xlabel('Time [s]')
ylabel('Pressure [bar]')
legend(fname)
%title(fname + " (" + cycles + " cycles)");

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


function [V, V_Theta] = volumeCycle(Ca)
R       = 8.5;
Vt      = 196 *10^(-6);                              % [m^2]
Vc      = Vt/R;                                 % [m^2]
Vd      = Vt-Vc;
phi     = 0;            %(360-153.4)/2/pi;

V=-Vd/2*cos(Ca*(2*pi/360))+Vc+Vd/2;

r = 0.030;
l = 0.085;
V_c = Vc
B = 0.06;
x = r*cos(Ca/360*2*pi) + sqrt(l^2 - r^2*(sin(Ca/360*2*pi))^2);

d_Theta = l + r - x;
V_Theta = pi*(B/2)^2*d_Theta +V_c;

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
