
function [V, V_Theta] = volumeCycle(Ca, Vt, R, phi)
Vc      = Vt/R;        
Vd      = Vt-Vc;

V=-Vd/2*cos(Ca*(2*pi/360))+Vc+Vd/2;

r = 0.030;
l = 0.085;
V_c = Vc
B = 0.06;
x = r*cos(Ca/360*2*pi) + sqrt(l^2 - r^2*(sin(Ca/360*2*pi))^2);

d_Theta = l + r - x;
V_Theta = pi*(B/2)^2*d_Theta +V_c;

end