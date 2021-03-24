clear all


data= readtable('E0_fl.xlsx');
pressure= data(:,2);
volume= data(:,1);

p= table2array(pressure);
v= table2array(volume);

figure(3)
hold on
plot(v,p)



xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')