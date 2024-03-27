function y = dopingFunction(x)

sinFunc1 = @(x) 249e15 * sin(20*pi*(x-0.075)) + 251e15;
sinFunc2 = @(x) 249e15 * sin(20*pi*(x-0.475)) + 251e15;


dopingDensity = @(x) (x>=0 & x<=0.1)*5e17 ...
    + (x>0.1 & x<0.15) .* sinFunc1(x) ...
    + (x>=0.15 & x<=0.45) * 2e15 ...
    + (x>0.45 & x<0.5) .* sinFunc2(x) ...
    + (x>=0.5 & x<=0.6) * 5e17;


y = dopingDensity(x);
end
