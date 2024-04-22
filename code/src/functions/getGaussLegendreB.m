function B = getGaussLegendreB(f,a,b,n)
if nargin ==3
    [~,B] = gaussLegendre(f,a,b);
elseif nargin == 4 
    [~,B] = gaussLegendre(f,a,b,n);
end
end