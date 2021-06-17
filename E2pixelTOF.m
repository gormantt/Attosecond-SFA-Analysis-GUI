function y = E2pixelTOF(x)
y0   = 325.59;
A1   = 5893.1;
t1   = 1.55726;
A2   = 4090.6;
t2   = 5.94996;
A3   = 3116.6;
t3   = 31.41142;
y = y0 + A1*exp(-x/t1) + A2*exp(-x/t2) + A3*exp(-x/t3);
end