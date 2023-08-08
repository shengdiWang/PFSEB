%
function snCap = snowThermalCap(roSnow)
% snow thermal capacity from HTESSEL

roIce = 920;
cIce  = 2.05;

snCap = roSnow * (cIce/roIce);       

end
