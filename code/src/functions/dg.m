function y = dg(x)
f = @(x) exp(-1./x) .*(x>0);
df = @(x) (x>0) .* exp(-1./x).*(1./x).^2 + (x==0) .* 0;
    y = zeros(size(x));
    for i = 1:length(x)
        if x(i) > 0 && x(i) < 1
            y(i) = (df(x(i)) .* f(1-x(i)) - f(x(i)) .* df(1-x(i))) ./ ((f(x(i)) + f(1-x(i))).^2);
        else
            y(i) = 0;
        end
    end
end
