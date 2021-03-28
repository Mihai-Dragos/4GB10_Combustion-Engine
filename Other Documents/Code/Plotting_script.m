clear all
%%% PLots same fuel different load conditions
% 
% data= readtable('E15_nl.xlsx');
% pressure= data(:,2);
% volume= data(:,1);
% 
% p= table2array(pressure);
% v= table2array(volume);
% 
% data2= readtable('E15_hl.xlsx');
% pressure2= data2(:,2);
% volume2= data2(:,1);
% 
% p2= table2array(pressure2);
% v2= table2array(volume2);
% 
% data3= readtable('E15_fl.xlsx');
% pressure3= data3(:,2);
% volume3= data3(:,1);
% 
% p3= table2array(pressure3);
% v3= table2array(volume3);
% 
% figure(3)
% grid on;
% grid minor;
% hold on
% plot(v,p)
% hold on
% plot(v2,p2)
% hold on
% plot(v3,p3)
% title("Plot of the pV Diagram for E_{15}, under specific load conditions ")
% legend("no load","half load","full load")
% 
% 
% xlabel('Volume [cm^3]')
% ylabel('Pressure [bar]')


%%
%%Plots same load different fuel type 

% data= readtable('E0_fl.xlsx');
% pressure= data(:,2);
% volume= data(:,1);
% 
% p= table2array(pressure);
% v= table2array(volume);
% 
% data2= readtable('E5_fl.xlsx');
% pressure2= data2(:,2);
% volume2= data2(:,1);
% 
% p2= table2array(pressure2);
% v2= table2array(volume2);
% 
% data3= readtable('E10_fl.xlsx');
% pressure3= data3(:,2);
% volume3= data3(:,1);
% 
% p3= table2array(pressure3);
% v3= table2array(volume3);
% 
% data4= readtable('E15_fl.xlsx');
% pressure4= data4(:,2);
% volume4= data4(:,1);
% 
% p4= table2array(pressure4);
% v4= table2array(volume4);
% 
% figure(3)
% grid on;
% grid minor;
% hold on
% plot(v,p)
% hold on
% plot(v2,p2)
% hold on
% plot(v3,p3)
% hold on
% plot(v4,p4)
% title("Plot of the pV Diagram for load condition: full load, under specific fuel types")
% legend("E0","E5","E10", "E15")
% 
% 
% xlabel('Volume [cm^3]')
% ylabel('Pressure [bar]')


%%
%%Plots Error Analysis All

data= readtable('E0_nl.xlsx');
pressure= data(:,2);
volume= data(:,1);
upper= data(:,3);
lower= data(:,4);
min= data(:,5);
max= data(:,6);

p= table2array(pressure);
v= table2array(volume);
u= table2array(upper);
l= table2array(lower);
mi= table2array(min);
ma= table2array(max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data2= readtable('E0_hl.xlsx');
pressure2= data2(:,2);
volume2= data2(:,1);
upper2= data2(:,3);
lower2= data2(:,4);
min2= data2(:,5);
max2= data2(:,6);

p2= table2array(pressure2);
v2= table2array(volume2);
u2= table2array(upper2);
l2= table2array(lower2);
mi2= table2array(min2);
ma2= table2array(max2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data3= readtable('E0_fl.xlsx');
pressure3= data3(:,2);
volume3= data3(:,1);
upper3= data3(:,3);
lower3= data3(:,4);
min3= data3(:,5);
max3= data3(:,6);

p3= table2array(pressure3);
v3= table2array(volume3);
u3= table2array(upper3);
l3= table2array(lower3);
mi3= table2array(min3);
ma3= table2array(max3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data4= readtable('E5_nl.xlsx');
pressure4= data4(:,2);
volume4= data4(:,1);
upper4= data4(:,3);
lower4= data4(:,4);
min4= data4(:,5);
max4= data4(:,6);

p4= table2array(pressure4);
v4= table2array(volume4);
u4= table2array(upper4);
l4= table2array(lower4);
mi4= table2array(min4);
ma4= table2array(max4);

data5= readtable('E5_hl.xlsx');
pressure5= data5(:,2);
volume5= data5(:,1);
upper5= data5(:,3);
lower5= data5(:,4);
min5= data5(:,5);
max5= data5(:,6);

p5= table2array(pressure5);
v5= table2array(volume5);
u5= table2array(upper5);
l5= table2array(lower5);
mi5= table2array(min5);
ma5= table2array(max5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data6= readtable('E5_fl.xlsx');
pressure6= data6(:,2);
volume6= data6(:,1);
upper6= data6(:,3);
lower6= data6(:,4);
min6= data6(:,5);
max6= data6(:,6);

p6= table2array(pressure6);
v6= table2array(volume6);
u6= table2array(upper6);
l6= table2array(lower6);
mi6= table2array(min6);
ma6= table2array(max6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data7= readtable('E10_nl.xlsx');
pressure7= data7(:,2);
volume7= data7(:,1);
upper7= data7(:,3);
lower7= data7(:,4);
min7= data7(:,5);
max7= data7(:,6);

p7= table2array(pressure7);
v7= table2array(volume7);
u7= table2array(upper7);
l7= table2array(lower7);
mi7= table2array(min7);
ma7= table2array(max7);

data8= readtable('E10_hl.xlsx');
pressure8= data8(:,2);
volume8= data8(:,1);
upper8= data8(:,3);
lower8= data8(:,4);
min8= data8(:,5);
max8= data8(:,6);

p8= table2array(pressure8);
v8= table2array(volume8);
u8= table2array(upper8);
l8= table2array(lower8);
mi8= table2array(min8);
ma8= table2array(max8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data9= readtable('E10_fl.xlsx');
pressure9= data9(:,2);
volume9= data9(:,1);
upper9= data9(:,3);
lower9= data9(:,4);
min9= data9(:,5);
max9= data9(:,6);


p9= table2array(pressure9);
v9= table2array(volume9);
u9= table2array(upper9);
l9= table2array(lower9);
mi9= table2array(min9);
ma9= table2array(max9);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data11= readtable('E15_nl.xlsx');
pressure11= data11(:,2);
volume11= data11(:,1);
upper11= data11(:,3);
lower11= data11(:,4);
min11= data11(:,5);
max11= data11(:,6);

p11= table2array(pressure11);
v11= table2array(volume11);
u11= table2array(upper11);
l11= table2array(lower11);
mi11= table2array(min11);
ma11= table2array(max11);

data12= readtable('E15_hl.xlsx');
pressure12= data12(:,2);
volume12= data12(:,1);
upper12= data12(:,3);
lower12= data12(:,4);
min12= data12(:,5);
max12= data12(:,6);

p12= table2array(pressure12);
v12= table2array(volume12);
u12= table2array(upper12);
l12= table2array(lower12);
mi12= table2array(min12);
ma12= table2array(max12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data13= readtable('E15_fl.xlsx');
pressure13= data13(:,2);
volume13= data13(:,1);
upper13= data13(:,3);
lower13= data13(:,4);
min13= data13(:,5);
max13= data13(:,6);

p13= table2array(pressure13);
v13= table2array(volume13);
u13= table2array(upper13);
l13= table2array(lower13);
mi13= table2array(min13);
ma13= table2array(max13);





t = tiledlayout(4,3); % Requires R2019b or later
t.Padding = 'compact';
t.TileSpacing = 'compact';


nexttile
grid on;
grid minor;
hold on
plot(v,p)
hold on
plot(v,u)
hold on
plot(v,l)
hold on
plot(v,mi)
hold on
plot(v,ma)
title("pV Diagram, error analysis E_{0}, no load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')


nexttile
grid on;
grid minor;
hold on
plot(v2,p2)
hold on
plot(v2,u2)
hold on
plot(v2,l2)
hold on
plot(v2,mi2)
hold on
plot(v2,ma2)
title("pV Diagram, error analysis E_{0}, half load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')


nexttile
grid on;
grid minor;
hold on
plot(v3,p3)
hold on
plot(v3,u3)
hold on
plot(v3,l3)
hold on
plot(v3,mi3)
hold on
plot(v3,ma3)
ylim([0 30])
title("pV Diagram, error analysis E_{0}, full load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')


nexttile
grid on;
grid minor;
hold on
plot(v4,p4)
hold on
plot(v4,u4)
hold on
plot(v4,l4)
hold on
plot(v4,mi4)
hold on
plot(v4,ma4)
title("pV Diagram, error analysis E_{5}, no load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')


nexttile
grid on;
grid minor;
hold on
plot(v5,p5)
hold on
plot(v5,u5)
hold on
plot(v5,l5)
hold on
plot(v5,mi5)
hold on
plot(v5,ma5)
title("pV Diagram, error analysis E_{5}, half load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')


nexttile
grid on;
grid minor;
hold on
plot(v6,p6)
hold on
plot(v6,u6)
hold on
plot(v6,l6)
hold on
plot(v6,mi6)
hold on
plot(v6,ma6)
title("pV Diagram, error analysis E_{5}, full load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')


nexttile
grid on;
grid minor;
hold on
plot(v7,p7)
hold on
plot(v7,u7)
hold on
plot(v7,l7)
hold on
plot(v7,mi7)
hold on
plot(v7,ma7)
title("pV Diagram, error analysis E_{10}, no load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')


nexttile
grid on;
grid minor;
hold on
plot(v8,p8)
hold on
plot(v8,u8)
hold on
plot(v8,l8)
hold on
plot(v8,mi8)
hold on
plot(v8,ma8)
title("pV Diagram, error analysis E_{10}, half load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')


nexttile
grid on;
grid minor;
hold on
plot(v9,p9)
hold on
plot(v9,u9)
hold on
plot(v9,l9)
hold on
plot(v9,mi9)
hold on
plot(v9,ma9)
title("pV Diagram, error analysis E_{10}, full load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')


nexttile
grid on;
grid minor;
hold on
plot(v11,p11)
hold on
plot(v11,u11)
hold on
plot(v11,l11)
hold on
plot(v11,mi11)
hold on
plot(v11,ma11)

title("pV Diagram, error analysis E_{15}, no load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')


nexttile
grid on;
grid minor;
hold on
plot(v12,p12)
hold on
plot(v12,u12)
hold on
plot(v12,l12)
hold on
plot(v12,mi12)
hold on
plot(v12,ma12)

title("pV Diagram, error analysis E_{15}, half load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')

nexttile
grid on;
grid minor;
hold on
plot(v13,p13)
hold on
plot(v13,u13)
hold on
plot(v13,l13)
hold on
plot(v13,mi13)
hold on
plot(v13,ma13)

title("pV Diagram, error analysis E_{15}, full load")
legend("Average plot", "Up appx. 2*sd", "Low appx. 2*sd","Minimum value", "Maximum value")
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')



