function y = dopingFunction(x)


sinFunc1 = @(x) 249 * sin(20*pi*(x-0.075)) + 251;
sinFunc2 = @(x) 249 * sin(20*pi*(x-0.475)) + 251;


% dopingDensity = @(x) (x>=0 & x<=0.1)*500 ...
%     + (x>0.1 & x<0.15) .* transition1(x) ...
%     + (x>=0.15 & x<=0.45) * 2 ...
%     + (x>0.45 & x<0.5) .* transition2(x) ...
%     + (x>=0.5 & x<=0.6) * 500;
% period
x = mod(x, 0.6);
% \mu m
y = dopingDensity(x) * 10e2;
end
function y = dopingDensity(x)
f = @(x) exp(-1./x) .*(x>0);
g = @(x) f(x) ./ (f(x) + f(1-x));

transition1 = @(x) 498 * g((0.15-x)/(0.15-0.1)) +2;
transition2 = @(x) 498 * g((x-0.45)/(0.5-0.45)) +2;
[m, n] = size(x);
y = zeros(m,n);
y = y(:);
for i=1:m*n
    if (x(i)>=0 && x(i)<=0.1)
        y(i) = 500;
    elseif x(i) > 0.1 && x(i) < 0.15
        y(i) = transition1(x(i));
    elseif x(i) >=0.15 && x(i) <=0.45
        y(i) = 2;
    elseif x(i)>0.45 && x(i) <0.5
        y(i) = transition2(x(i));
    elseif x(i)>=0.5 && x(i)<= 0.6
        y(i) = 500;
    end
end
y = reshape(y,m,n);
end
