clear all;clc
%% add general to matlab path 
addpath('General');
%%
DataDir='C:\Users\20190762\Documents\Universiteit\2020-2021 (Q3) Combustion Engine (4GB10)\4GB10_matlab_and_trainingdata\matlab\Data\No_fuel';ColumnOrder={'time','Sensor','Encoder'};
% DataDir='../Data/Gasoline';cOrder={'time'','Encoder','Sensor};
col = lines(3);
%% Loading all measurments in DataDir
% This script can be used to determine the position of the double tooth
% position using max pressure of a No-fuel trial. 
% A lot of the data collected in first two attempts is poor, and useless. 
% Small error analysis is provided at the end

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
    V      = Data.Volt;
    RevEnd = Data.RevEnds;
    
    if Data.NRevs >= 2;
        [~, maxPressureID] = findpeaks(Data.Volt, 'MinPeakHeight', 0.5, 'MinPeakDistance', 3000)
        stepAngle(i) = 360/(abs(Data.RevEnds(2) - Data.RevEnds(3)));
        for ii = 1:2
           stepDiff=abs(maxPressureID(2) - Data.RevEnds(ii));
           angle(i,ii) =stepAngle(i)*stepDiff;
        end
    else
        continue 
    end
end

%% Error Anaylsis

meanAngle = mean(angle);
stdAngle = std(angle);

%%