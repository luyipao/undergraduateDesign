function [obj, E, func] = GK3(obj)
t = obj.CFL * obj.meshSize^2;
k0 = obj.coeffs;
x = linspace(0,0.6,1000);
F = L(obj);
k1 = k0 + t * F;
obj.coeffs = k1;
F = L(obj);
k2 = 3/4 * k0 + 1/4 * k1 + 1/4 * t * F;
obj.coeffs = k2;
F = L(obj);
obj.coeffs = 1/3 * k0 + 2/3 * k2 + 2/3 * t * F;
[~,E,func] = L(obj);
end