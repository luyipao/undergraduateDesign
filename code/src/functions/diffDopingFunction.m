function y = diffDopingFunction(x)

% y = 20 * pi * 249 * (x>= 0.1 & x<= 0.15) .* cos(20*pi*(x-0.075)) ...
%     + 20 * pi * 249 * (x>= 0.45 & x<= 0.5) .* cos(20*pi*(x-0.475));
% y = y * 10e2;
% f = @(x) exp(-1./x) .*(x>0);
% g = @(x) (f(x) ./ (f(x) + f(1-x))) .* (x>0 & x<1);
% df = @(x) (x>0) .* exp(-1./x).*(1./x).^2 + (x==0) .* 0;
% dg = @(x) ((df(x) .* f(1-x) - f(x).*df(1-x)) ./ (f(x) + f(1-x)).^2) .* (x>0 & x<1);
%

transition1 = @(x)  -9960 * dg((0.15-x)/(0.15-0.1));
transition2 = @(x) 9960 *dg((x-0.45)/(0.5-0.45));

dopingDensity = @(x)  (x>0.1 & x<0.15) .* transition1(x) ...
    + (x>0.45 & x<0.5) .* transition2(x);
% period
x = mod(x, 0.6);
% \mu m
y = dopingDensity(x) * 10e2;
end

function y = dg(x)
f = @(x) exp(-1./x) .*(x>0);
df = @(x) (x>0) .* exp(-1./x).*(1./x).^2 + (x==0) .* 0;
y = zeros(size(x));
for i = 1:length(x)
    if x(i) >0 && x(i) < 1
        y(i) = (df(x(i)) .* f(1-x(i)) + f(x(i)) .* df(1-x(i))) ./ ((f(x(i)) + f(1-x(i))).^2);
    else
        y(i) = 0;
    end
end
end