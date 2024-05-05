clear;
f = @(x) dopingFunction(x);
a = 0; b = 0.6; n = 2;iterNum = 5;

L2Error = [];
mesh = 25;
N = mesh;
for i = 1:iterNum
     X = linspace(a,b,1000);
     fPiecewiseProj = piecewiseL2Projection(f,n,a,b,mesh * 2^(i-1));
     L2Error = [L2Error sqrt( sum(abs(fPiecewiseProj(X)' - f(X)).^2) )];
end

mesh = 12 .* 2.^(0:iterNum-1);
Order = [NaN log2(L2Error(1:end-1) ./ L2Error(2:end))];  
% table
T = table(mesh', L2Error', Order', 'VariableNames', {'mesh', 'L2Error', 'Order'});
T.L2Error = num2str(T.L2Error, '%.2e');

X = linspace(a,b,1000);
fPiecewiseProj = piecewiseL2Projection(f,n,a,b,96);
plot(X,f(X),'-',X,fPiecewiseProj(X),'--');
legend('exact','piece wise projection');

% 
