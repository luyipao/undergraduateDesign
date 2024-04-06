addpath('C:\Users\12192\OneDrive - zju.edu.cn\documents\homeworks\undergraduateDesign\code\src\functions');

a = 0;
b = 0.6;
n = 3;
N = 100;

mesh = Mesh(a,b,N,n);
mesh = mesh.ENOreconstruction(@(x) dopingFunction(x));
mesh = mesh.auxiliaryDDModelDGFunction(@(x) dopingFunction(x));
mesh.draw();
% [dopingProj,dopingProjCoeffVec] = piecewiseL2Projection(@(x) dopingFunction(x),n,a,b,N);
% [dopingProj] = getElectronConcentration(dopingProjCoeffVec, mesh, n);
% [auxq,auxiliaryVarCoeffVec] = auxiliaryDDModelDGFunction(mesh,@(x) dopingFunction(x),n);


% epsilon = 0.001;
% CFL =1/(2*n+1);
% [C, Q, X, T] = DDModelDGFunction(dopingProjCoeffVec,@(x) sqrt(THETA * RELAXATION_PARAMETER) * diffDopingFunction(x), mesh, n, epsilon, CFL);
% 

% %