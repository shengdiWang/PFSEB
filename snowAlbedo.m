%
function snalbedo = snowAlbedo(roSnow)
% snow albedo from Eq. 31

snalbedo = 1 - 0.247 * sqrt(0.16 + 110 * (roSnow/1000)^4);

end

function snTC = snowThermalCon(roSnow)
% snow thermal conductivity from HTESSEL

roIce = 920;
kIce  = 2.29;

snTC = kIce * (roSnow / roIce)^1.88; % snow conductivity

end