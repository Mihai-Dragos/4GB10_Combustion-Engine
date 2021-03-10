close all;
dCa(1) = 1
Ca(1) = 0
[V1(1), V2(1)] = volumeCycle(Ca(1));
for i = 2:361
    Ca(i) = Ca(i-1) +dCa;
    [V1(i), V2(i)] = volumeCycle(Ca(i));
    dV(i) = V2(i) - V1(i);
end
figure()
hold on
plot(dV)
plot(V1)
plot(V2)
legend("Volume Differnece", "Our Formula", "Their Formula")
xlabel("Crrank angle [deg]")
ylabel("Volume in cylinder [m^3]")
grid on
grid minor

