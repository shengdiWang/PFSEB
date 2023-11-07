%
function dLWout = derivativeLongWaveOut(emiV, Ts0)
% snow thermal capacity from HTESSEL

SIGMA = 5.670374E-8;  % Stefan-Boltzmann constant [W m-2 K-4]
dLWout = -4.0 * emiV * SIGMA * Ts0^3;

end