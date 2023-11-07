%
function snTC = snowThermalCon(roSnow)
% snow thermal conductivity from HTESSEL

roIce = 920;
kIce  = 2.29;

snTC = kIce * (roSnow / roIce)^1.88; % snow conductivity

end