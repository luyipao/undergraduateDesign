
addpath('.\functions');
addpath('..\src\functions')
setParameters

mesh = linspace(a,b,N+1);
meshSize = (b-a)/N;

[dopingProj,dopingProjCoeffVec] = piecewiseL2Projection(@(x) dopingFunction(x),n,a,b,N);

[auxq,auxiliaryVarCoeffVec] = auxiliaryDDModelDGFunction(mesh,dopingProj,n);


epsilon = 1;
CFL = 0.2;
[C, Q, X, T] = DDModelDGFunction(dopingProjCoeffVec,auxq, mesh, n, epsilon, CFL);
