

addpath('C:\Users\12192\OneDriv e - zju.edu.cn\documents\homeworks\undergraduateDesign\code\src\functions');

setParameters
a = 0;
b = 0.6;
n = 3;
N = 96;

mesh = linspace(a,b,N+1);
meshSize = (b-a)/N;

[dopingProj,dopingProjCoeffVec] = piecewiseL2Projection(@(x) dopingFunction(x),n,a,b,N);

[auxq,auxiliaryVarCoeffVec] = auxiliaryDDModelDGFunction(mesh,dopingProj,n);


epsilon = 0.001;
CFL =1/(2*n+1);
[C, Q, X, T] = DDModelDGFunction(dopingProjCoeffVec,@(x) sqrt(THETA * RELAXATION_PARAMETER) * diffDopingFunction(x), mesh, n, epsilon, CFL);


%  