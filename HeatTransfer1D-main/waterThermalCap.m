%
function waterCap = waterThermalCap(T)
% water thermal capacity from HTESSEL as Eq. 22

waterCap = 4.20843 + 1.11362E-1 * T + 5.12142E-3 * T^2 + 9.3482E-5 * T^3; 

end