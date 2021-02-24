%%
%%For full_load dataset


<<<<<<< Updated upstream
full_table = dlmread("C:\Users\alexa\OneDrive\Desktop\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\full_load.txt");
full_time = full_table(:,1);
full_voltage = full_table(:,2);
full_pulse = full_table(:,3);
%full_volume = full_table(:,4);
=======
full_table = dlmread("C:\Users\20161501\OneDrive - TU Eindhoven\Year3\Q3\4GB10 DBL Combustion engine\GITKRAKEN\4GB10_Combustion-Engine\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\full_load.txt");
full_time = full_table(:,1);
full_voltage = full_table(:,2);
full_pulse = full_table(:,3);


>>>>>>> Stashed changes



V_s = 5; %V, constant

V_full = full_voltage * V_s;
%2nd voltage 
%3rd pulse

%1 complete cycle is 0.04s = 1/25

P_full = ((V_full/V_s-0.115)/0.0385); 

figure()
hold on
title("full")
plot(full.t,P_full)
xlim([0.035 0.075])
hold off

%%
%%For half_load dataset
<<<<<<< Updated upstream

half_table = dlmread("C:\Users\alexa\OneDrive\Desktop\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\half_load.txt");

=======
half_table = dlmread("C:\Users\20161501\OneDrive - TU Eindhoven\Year3\Q3\4GB10 DBL Combustion engine\GITKRAKEN\4GB10_Combustion-Engine\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\half_load.txt");
>>>>>>> Stashed changes
half_time = half_table(:,1);
half_voltage = half_table(:,2);
half_pulse = half_table(:,3);



V_s = 5; %V, constant

V_half = half_voltage * V_s;
%2nd voltage 
%3rd pulse

%1 complete cycle is 0.04s = 1/25

P_half = ((V_half/V_s-0.115)/0.0385); 

figure()
hold on
title("half")
plot(full.t,P_half)
xlim([0 0.04])
hold off




%%
%%For no_load dataset
no_table = dlmread("C:\Users\20161501\OneDrive - TU Eindhoven\Year3\Q3\4GB10 DBL Combustion engine\GITKRAKEN\4GB10_Combustion-Engine\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\no_load.txt");

<<<<<<< Updated upstream

no_table = dlmread("C:\Users\alexa\OneDrive\Desktop\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\no_load.txt");


=======
>>>>>>> Stashed changes
no_time = no_table(:,1);
no_voltage = no_table(:,2);
no_pulse = no_table(:,3);



V_s = 5; %V, constant


V_no = no_voltage * V_s;%2nd voltage 
%3rd pulse

%1 complete cycle is 0.04s = 1/25

P_no = ((V_no/V_s-0.115)/0.0385); 
figure()
hold on
title("no")
plot(full.t,P_no)
xlim([0 0.04])
hold off


<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes
%% Formula for angle and volume
%t_theta = 

theta_0 = 180; %crank angle at te begin of a cycle [degree]
v_theta= 18000; %[degree/second]
k=1; 
for i=0:full_time
    theta = v_theta*full_time(i) + theta_0;
if theta > 360
    theta = theta - 360*k;
    k=k+1;
end
end

% d_tetha = l+r-r*cos(tetha)-sqrt(l^2-r^2*sin*(tetha)^2);

%V_tetha = pi*(B/2)^2*d_tetha+V_c;

%% Pulse full load

% figure()
% hold on
% plot(full_time, full_pulse)
% hold off

%% Pulse

figure()
hold on
plot(full_time, full_pulse)
hold off
<<<<<<< Updated upstream
=======

>>>>>>> Stashed changes

figure()
hold on
plot(full_time, full_pulse)

title("Full load pulse sensor");
xlim([0 0.04])
hold off
%% Pulse graphs half load
% figure()
% hold on
% plot(half_time, half_pulse)
% hold off

figure()
hold on
plot(half_time, half_pulse)
title("Half load pulse sensor");
xlim([0 0.04])
hold off

%% Pulse graphs no load
% figure()
% hold on
% plot(no_time, no_pulse)
% hold off

xlim([0 0.04])
hold off
%%
figure()
hold on
plot(full_time, full_pulse)
hold off

figure()
hold on
plot(full_time, full_pulse)
xlim([0 0.04])
hold off

%%
figure()
hold on
plot(no_time, no_pulse)
hold off


figure()
hold on
plot(no_time, no_pulse)

title("No load pulse sensor");
xlim([0 0.04])
hold off

%% Pressure graphs no load
% figure()
% hold on
% plot(no_time, P_no)
% hold off

xlim([0 0.04])
hold off

%%
figure()
hold on
plot(no_time, P_no)
hold off


figure()
hold on
plot(no_time, P_no)
xlim([0 0.04])

title("No load pressure sensor");
hold off

%% Pressure graphs half load
% figure()
% hold on
% plot(half_time, P_half)
% hold off

figure()
hold on
x=1;
plot(half_time, P_half)
plot(x)
xlim([0 0.04])
title("Half load pressure sensor");
hold off

%% Pressure graphs full load

% figure()
% hold on
% plot(full_time, P_full)
% hold off

figure()
hold on
plot(full_time, P_full)
xlim([0 0.04])
title("Full load pressure sensor");
<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes
hold off

%%
figure()
hold on
plot(half_time, P_half)
hold off

figure()
hold on
plot(half_time, P_half)
xlim([0 0.04])
hold off


%%
figure()
hold on
plot(full_pulse,P_full)
xlim([0 0.04])
<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes
hold off


%%
<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes
% figure()
% hold on
% plot(full_pulse,P_full)
% xlim([0 0.04])
% hold off
% 
% 
%%
% figure()
% hold on
% plot(full.Vol,P_full)
<<<<<<< Updated upstream
%  hold off

=======
% 
% hold off
>>>>>>> Stashed changes
figure()
hold on
plot(full.Vol,P_full)

hold off
<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes


