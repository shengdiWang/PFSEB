%
function thaov = vthao(sand)

slopev     = -0.0095;
interceptv = 1.54;
thaov      = slopev * sand + interceptv;
thaov      = 10^thaov / 1000; %[m]

end