%
function avp = atmosphericVaporPressure(tDew)
% Atmospheric vapor pressure as Eq. 3

%     avp = 10^(11.40 - 2353.0 / tDew);
    tDew = tDew - 273.15;
    avp = 611.2 * exp(17.62 * tDew/(243.12 + tDew));
end