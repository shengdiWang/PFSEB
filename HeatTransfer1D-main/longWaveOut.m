%
function LWout = longWaveOut(emiV, Ts0)
% Emitted longwave radiation as Eq. 4

SIGMA = 5.670374E-8;  % Stefan-Boltzmann constant [W m-2 K-4]
LWout = - emiV * SIGMA * Ts0 ^ 4;

end