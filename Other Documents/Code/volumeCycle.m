
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