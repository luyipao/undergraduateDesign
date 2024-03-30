function y = dopingFunction(x)

sinFunc1 = @(x) 249 * sin(20*pi*(x-0.075)) + 251;
sinFunc2 = @(x) 249 * sin(20*pi*(x-0.475)) + 251;


dopingDensity = @(x) (x>=0 & x<=0.1)*500 ...
    + (x>0.1 & x<0.15) .* sinFunc1(x) ...
    + (x>=0.15 & x<=0.45) * 2 ...
    + (x>0.45 & x<0.5) .* sinFunc2(x) ...
    + (x>=0.5 & x<=0.6) * 500;

% \mu m
y = dopingDensity(x) * 10e2;
end
