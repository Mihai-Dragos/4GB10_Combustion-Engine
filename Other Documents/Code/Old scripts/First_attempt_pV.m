Tv1 = linspace(273, 323, 10);
Tv2 = linspace(323, 373, 10);
Vv1 = linspace(150, 200, 10);
[T1m,V1m] = ndgrid(Tv1,Vv1);
[T2m,V1m] = ndgrid(Tv2,Vv1);
R = 0.2871;
m = 5;
% p1 = 7.8378;
% V1 = 50;
p = @(T,V) m*R*T./V;
figure
surf(T1m, V1m, p(T1m,V1m))
hold on
surf(T2m, V1m, p(T2m,V1m))
hold off
grid on
xlabel('T')
ylabel('V')
zlabel('p')
view(90,0)
figure
plot(Vv1, p(Tv1,Vv1), 'LineWidth',2)
hold on
plot(Vv1, p(Tv2,Vv1), 'LineWidth',2)
plot(Vv1(1)*[1 1], p([Tv1(1) Tv1(end)],[1 1]*Vv1(1)), '-g', 'LineWidth',2)
plot(Vv1(end)*[1 1], p([Tv2(1) Tv2(end)],[1 1]*Vv1(end)), '-g', 'LineWidth',2)
hold off
grid
xlabel('V')
ylabel('p')
axis([125  225    2.0  3.3])
legend('T(273 - 323)', 'T(323 - 373)')