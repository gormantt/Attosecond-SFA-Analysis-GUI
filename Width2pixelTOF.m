function y = Width2pixelTOF(x)
y0   = 13;
A1   = 3.5e8;
xc   = -21;
P    = -4.45;
y = y0 + A1*(x-xc)^P;
end