clear;

% 确定待处理函数，区间及初始数量
f = @(x)  sin(x);
a = 0; b = 2*pi; n = 3;

L2Error = [];
mesh = 10;
N = mesh;
for i = 1:5
     X = linspace(a,b,1000);
     fPiecewiseProj = piecewiseL2Projection(f,n,a,b,mesh * 2^(i-1));
     L2Error = [L2Error sqrt(norm(fPiecewiseProj(X) - f(X))^2)];
end

mesh = 10 .* 2.^(0:4);
Order = [NaN log2(L2Error(1:end-1) ./ L2Error(2:end))];  
% table
T = table(mesh', L2Error', Order', 'VariableNames', {'mesh', 'L2Error', 'Order'});
T.L2Error = num2str(T.L2Error, '%.2e');
T