
addpath('scripts');
addpath('functions');
setParameters
DD = Model('input.json');
a = 0;
b = 0.6;
n = 2;
nSteps = zeros(2,5);
N = [100 200];
t = [1.6e-5];
M = cell(length(N), length(t));
for i = 1:length(N)
    for j = 1:length(t)
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, DD);
        mesh.t = t(j);
        [f, M{i,j}] = mesh.DDModelDGFunction(2);
    end
end
x = linspace(0,0.6,1000);
        
% figure(1)
% x = linspace(0,0.6,1000);
% electronConcentration = obj.getBasisPolys(obj.coeffs);
% plot(x,electronConcentration.solve(x));
% figure(2)
% % Generate a filename with N and t
% E = obj.getBasisPolys(obj.Ecoeffs);
% plot(x,E.solve(x));

% k = 4;
% Funcs = cell(k,1);
% L2Error = [];
%
% for i=1:k
%     x = linspace(0,0.6);
%     N = 25;
%     mesh = Mesh(a,b,N*2^(i-1),@(x) dopingFunction(x),n);
%     Funcs{i} = mesh.DDModelDGFunction(n);
%     %L2Error = [L2Error sqrt( sum(abs(Funcs{i}.solve(x)' - dopingFunction(x)).^2) )];
% end
% Order = [NaN log2(L2Error(1:end-1) ./ L2Error(2:end))];



%test L2 projection
% f = piecewiseL2Projection(@(x) dopingFunction(x) ,n,a,b,N);
% x = linspace(0,0.6,1000);
% scatter(x,dopingFunction(x)-f(x)');
% plot(x,dopingFunction(x), x,f(x),'--');
% legend('exact', 'projection');


