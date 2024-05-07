
addpath('scripts');
addpath('functions');
setParameters
a = 0;
b = 0.6;
n = 2;
nSteps = zeros(2,5);
N = [100 200];
t = [1.6e-5 4.2e-6];
nmodel = 3;
% f = cell(length(N), length(t), nmodel);
% E = cell(length(N), length(t), nmodel);
% M = cell(nmodel,1);
% TV = cell(length(N), length(t), nmodel);
% e = cell(length(N), length(t), nmodel);
DD = Model('mobilityEqConst.json');
for i = 1:length(N)
    for j = i
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, DD);
        mesh.t = t(j);
        [tempM{i,j}, TV{i,j,1}, e{i,j,1}] = mesh.DDModelDGFunction(2);
    end
end
M{1} = tempM;
clear tempM;

N = 100;

DD2 = Model('mobilityDpDopingFunction.json');
for i = 1:length(N)
    for j = i
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, DD2);
        mesh.t = t(j);
        [tempM{i,j}, TV{i,j,2}, e{i,j,2}] = mesh.DDModelDGFunction(2);
    end
end
M{2} = tempM;
clear tempM;

HF = Model('HFmobilityEqConst.json');
for i = 1:length(N)
    for j = i
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, HF);
        mesh.t = t(j);
        [tempM{i,j}, TV{i,j,3}, e{i,j,3}] = mesh.DDModelDGFunction(2);
    end
end
M{3} = tempM;
clear tempM;

DD = Model('mobilityEqConst.json');
for i = 1:length(N)
    for j = i
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, DD);
        mesh.t = t(j);
        [tempM{i,j}, TV{i,j,4}, e{i,j,4}] = mesh.DDModelDGFunction(2);
    end
end
M{4} = tempM;
clear tempM;

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


