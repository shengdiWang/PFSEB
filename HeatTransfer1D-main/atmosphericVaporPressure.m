%
function avp = atmosphericVaporPressure(tDew)
% Atmospheric vapor pressure as Eq. 3

%     avp = 10^(11.40 - 2353.0 / tDew);
    avp = 611 * exp(17.62 * (tDew - 273.15) / (243.12 + (tDew - 273.15)));
end


