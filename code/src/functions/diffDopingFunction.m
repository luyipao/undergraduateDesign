function y = diffDopingFunction(x)

y = 20 * pi * 249 * (x>= 0.1 & x<= 0.15) .* cos(20*pi*(x-0.075)) ...
    + 20 * pi * 249 * (x>= 0.45 & x<= 0.5) .* cos(20*pi*(x-0.075));
y = y * 10e2
