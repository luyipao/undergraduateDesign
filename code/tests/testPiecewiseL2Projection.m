clear;
% 确定待处理函数，区间及初始数量
f = @(x) dopingFunction(x) / 10e17;
a = 0; b = 0.6; n = 3;iterNum = 6;

L2Error = [];
mesh = 12;
N = mesh;
for i = 1:iterNum
     X = linspace(a,b,1000);
     fPiecewiseProj = piecewiseL2Projection(f,n,a,b,mesh * 2^(i-1));
     L2Error = [L2Error sqrt(sum(abs(fPiecewiseProj(X) - f(X)).^2) )];
end

mesh = 12 .* 2.^(0:iterNum-1);
Order = [NaN log2(L2Error(1:end-1) ./ L2Error(2:end))];  
% table
T = table(mesh', L2Error', Order', 'VariableNames', {'mesh', 'L2Error', 'Order'});
T.L2Error = num2str(T.L2Error, '%.2e');
T

X = linspace(a,b,10000);
fPiecewiseProj = piecewiseL2Projection(f,n,a,b,96);
plot(X,f(X),'-',X,fPiecewiseProj(X),'--');
legend('exact','piece wise projection');