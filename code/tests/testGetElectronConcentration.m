%% testGetElectronConcentration
clear;
n = 3;
a = 0;
b = 0.6;
N = 12;


[dopingProj,dopingProjCoeffVec] = piecewiseL2Projection(@(x) dopingFunction(x),n,a,b,N);

mesh = linspace(a,b,N+1);

electronConcentration = getElectronConcentration(dopingProjCoeffVec,mesh,n);

X = linspace(a,b);

% draw
plot(X,dopingProj(X),'-',X,electronConcentration(X),'--');

legend('exact','vecGeneration');
