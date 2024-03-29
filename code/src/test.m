
addpath('C:\Users\12192\OneDrive - zju.edu.cn\documents\homeworks\undergraduateDesign\code\src\functions');

setParameters

mesh = linspace(a,b,N+1);
meshSize = (b-a)/N;

[dopingProj,dopingProjCoeffVec] = piecewiseL2Projection(@(x) dopingFunction(x),n,a,b,N);

[auxq,auxiliaryVarCoeffVec] = auxiliaryDDModelDGFunction(mesh,dopingProj,n);


epsilon = 0.001;
CFL = 0.2;
[C, Q, X, T] = DDModelDGFunction(dopingProjCoeffVec,auxq, mesh, n, epsilon, CFL);


% clear;
% ff = @(x) x;
% gg = @(x) x;
% f = @(x) arrayfun(@(x) ff(x), x);
% g = @(x) arrayfun(@(x) gg(x).^2, x);
% h = @(x) arrayfun(@(x) quadgk(@(y) f(y)+g(y), 0, x), x);
