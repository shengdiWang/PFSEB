%
function bv = vb(clay)

slopev     = 0.157;
interceptv = 3.10;
bv         = slopev * clay + interceptv;

end