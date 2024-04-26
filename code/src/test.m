
addpath('C:\Users\12192\OneDrive - zju.edu.cn\documents\homeworks\undergraduateDesign\code\src\functions');
setParameters
a = 0;
b = 0.6;
n = 3;
N = 100;
mesh = Mesh(a,b,N,@(x) dopingFunction(x),n);
mesh = mesh.DDModelDGFunction;

% k = 4;
% Funcs = cell(k,1);
% L2Error = [];
% 
% for i=1:k
%     x = linspace(0,0.6);
%     N = 25;
%     mesh = Mesh(a,b,N*2^(i-1),@(x) dopingFunction(x),n);
%     Funcs{i} = mesh.DDModelDGFunction;
%     L2Error = [L2Error sqrt( sum(abs(Funcs{i}.solve(x)' - dopingFunction(x)).^2) )];
% end
% Order = [NaN log2(L2Error(1:end-1) ./ L2Error(2:end))];  



%test L2 projection
f = piecewiseL2Projection(@(x) dopingFunction(x) ,n,a,b,N);
x = linspace(0,0.6,1000);
scatter(x,dopingFunction(x)-f(x)');
plot(x,dopingFunction(x), x,f(x),'--');
legend('exact', 'projection');