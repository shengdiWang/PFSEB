%
function wu = NiuYang(satW, t, b, taos)
% Niu & Yang:
%   Effects of Frozen Soil on Snowmelt Runoff and Soil Water Storage at a
%   Continental Scale
%   Eq.3
%
% The impact of new land surface physics on the GCM simulation of climate
% and climate sensitivity
% Eq. 30
%
% t [K]

t0 = 273.15;
lf = 334 * 1000; % [J/kg];
g  = 9.80665;

k  = lf / g / t;

a  = k * (t - t0) / taos;
wu = satW * (-1 * a)^(-1 / b);

end

