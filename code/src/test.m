
addpath('scripts');
addpath('functions');
setParameters
DD = Model('mobilityEqConst.json');
a = 0;
b = 0.6;
n = 2;
nSteps = zeros(2,5);
N = [100 200];
t = [1.6e-5 4.2e-6];
nmodel = 3;
f = cell(length(N), length(t), nmodel);
E = cell(length(N), length(t), nmodel);
M = cell(length(N), length(t), nmodel);
TV = cell(length(N), length(t), nmodel);
e = cell(length(N), length(t), nmodel);
for i = 1:length(N)
    for j = i
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, DD);
        mesh.t = t(j);
        [M{i,j}, TV{i,j}, e{i,j}] = mesh.DDModelDGFunction(2);
    end
end
x = linspace(0, 0.6);
[f{1,1,1}, E{1,1,1}] = M{1}{1,1}.getFuncs;
[f{2,2,1}, E{2,2,1}] = M{1}{2,2}.getFuncs;
%
figure(1)
plot(x, f{1,1,1}.solve(x), x, f{2,2,1}.solve(x), 'o', 'LineWidth', 1);
legend('N=100', 'N=200', 'Location', 'north', 'FontSize', 11);
xlabel('x');
ylabel('n');

gca.XMinorTick = 'on';
gca.YMinorTick = 'on';
figure(2)
plot(x, E{1,1,1}(x), x, f{2,2,1}(x), 'o', 'LineWidth', 1);
legend('N=100', 'N=200', 'Location', 'north', 'FontSize', 11);
xlabel('x');
ylabel('E');

gca.XMinorTick = 'on';
gca.YMinorTick = 'on';
% 
[f{1,1,2}, E{1,1,2}] = M{2}{1,1}.getFuncs;
[f{2,2,2}, E{2,2,2}] = M{2}{2,2}.getFuncs;

figure(3)
legend('\mu = 0.75', '\mu = \mu(n_d)', 'Location', 'north', 'FontSize', 11);
xlabel('x');
ylabel('n');
gca.XMinorTick = 'on';
gca.YMinorTick = 'on';

figure(4)
plot(x, E{1,1,1}(x), x, E{1,1,2}(x), 'o', 'LineWidth', 1);
legend('\mu = 0.75', '\mu = \mu(n_d)', 'Location', 'north', 'FontSize', 11);
xlabel('x');
ylabel('E');
gca.XMinorTick = 'on';
gca.YMinorTick = 'on';

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


