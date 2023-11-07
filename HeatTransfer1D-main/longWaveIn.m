%
function LWin = longWaveIn(EA, TA)
% Long-wave incoming radiation as Eq. 2

SIGMA = 5.670374E-8;  % Stefan-Boltzmann constant [W m-2 K-4]
LWin  = 1.08 * (1 - exp(-(0.01 * EA)^(TA / 2016))) * SIGMA * TA^4;

end
