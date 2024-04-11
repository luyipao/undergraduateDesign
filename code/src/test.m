% f1 = @(x) 500+ 0 * x;
% f2 = @(x) 2 + 0*x;
% f = @(x) exp(-1./x) .*(x>0);
% g = @(x) f(x) ./ (f(x) + f(1-x));
% a = 0.1;
% b = 0.15;
% c = 0.2;
% d = 0.25;
% X = linspace(a,d,1000);
% 
% plot(X,g( (X-a)/(b-a)) .* g((d-X)/(d-c)) )
% hold off
addpath('C:\Users\12192\OneDrive - zju.edu.cn\documents\homeworks\undergraduateDesign\code\src\functions');
setParameters
a = 0;
b = 0.6;
n = 3;
N = 100;

mesh = Mesh(a,b,N,@(x) dopingFunction(x),n);
mesh = mesh.DDModelDGFunction();

% [dopingProj,dopingProjCoeffVec] = piecewiseL2Projection(@(x) dopingFunction(x),n,a,b,N);
% [dopingProj] = getElectronConcentration(dopingProjCoeffVec, mesh, n);
% [auxq,auxiliaryVarCoeffVec] = auxiliaryDDModelDGFunction(mesh,@(x) dopingFunction(x),n);


% epsilon = 0.001;
% CFL =1/(2*n+1);
% [C, Q, X, T] = DDModelDGFunction(dopingProjCoeffVec,@(x) sqrt(THETA * RELAXATION_PARAMETER) * diffDopingFunction(x), mesh, n, epsilon, CFL);
% 

% %