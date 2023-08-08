%
function Lstar   = Monin_obokhov(roAir, Ustar, TA, QH, QE)

CAir = 1005.7;
L_w  = 1000 * (2500.8 - 2.36 * (TA - 273.15));   % [J/kg] latent heat of evaporation of water
VONK = 0.4;
GRAVIT = 9.81;


Lstar  = - roAir * CAir * Ustar^3 * TA / (VONK * GRAVIT * (QH + 0.61*CAir*TA*QE/L_w));
% fprintf('%6.2f\n', Lstar);

% Lstar = (abs(Lstar) < 1e-6) * Lstar / abs(Lstar) * 1e-6 + (abs(Lstar) >= 1e-6) * Lstar;  % lower limit for Lstar
% Lstar = (abs(Lstar) > 1e6) *  Lstar / abs(Lstar) * 1e6  + (abs(Lstar) <= 1e6)  * Lstar;     % upper limit for Lstar

end