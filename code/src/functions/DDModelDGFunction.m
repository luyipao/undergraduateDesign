%output n^h coeff matrix: C; n^h Value matrix electronConcentrationValue;
% space mesh: mesh; Time mesh: T.
function [C, electronConcentrationValue, mesh, T] = DDModelDGFunction(electronConcentrationCoeffVec,auxq, mesh, n, epsilon, CFL)
N = length(mesh) - 1;
meshSize = max(abs(mesh(2:end)-mesh(1:end-1)));

[~, ~, PiPjIntegral] = getLegendreBasisInfo(n,meshSize);
A = PiPjIntegral;

t = CFL * meshSize;

C = [electronConcentrationCoeffVec GK3(electronConcentrationCoeffVec,A,t,mesh,n)] ;
T = [0];
nLast = getElectronConcentration(C(:,end-1),mesh,n);
nNow = getElectronConcentration(C(:,end),mesh,n);
electronConcentrationValue = [nLast(mesh)' nNow(mesh)'];

while(norm(electronConcentrationValue(:,end-1) - electronConcentrationValue(:,end)) > epsilon)
    T = [T T(end)+t];
    C = [C GK3(electronConcentrationCoeffVec,A,t,mesh,n)] ;
    nNow = getElectronConcentration(C(:,end),mesh,n);
    electronConcentrationValue = [electronConcentrationValue nNow(mesh)'];
end
end





function C = GK3(C,A,t,mesh,n)
global F auxq E electronConcentration
electronConcentration = getElectronConcentration(C, mesh, n);
auxq = auxiliaryDDModelDGFunction(mesh, electronConcentration, n);
E = getElectricField(electronConcentration);
F = getF(mesh, auxq, electronConcentration, E, n);

k1 = C + A\F;

electronConcentration = getElectronConcentration(k1, mesh, n);
auxq = auxiliaryDDModelDGFunction(mesh, electronConcentration, n);
E = getElectricField(electronConcentration);
F = getF(mesh, auxq, electronConcentration, E, n);

k2 = 3/4 * C + 1/4 * k1 + A\F;

electronConcentration = getElectronConcentration(k2, mesh, n);
auxq = auxiliaryDDModelDGFunction(mesh, electronConcentration, n);
E = getElectricField(electronConcentration);
F = getF(mesh, auxq, electronConcentration, E, n);

C = 1/3 * C + 2/3 * k2 + 2/3 * t * A\F;
end

%%
function electricField = getElectricField(electronConcentration, C, mesh)
global ELECTRON_CHARGE DIELECTRIC_PERMITTIVITY
coeff = ELECTRON_CHARGE / DIELECTRIC_PERMITTIVITY;
electronConcentrationFunc = @(x) arrayfun(@(x) electronConcentration(mod(x,mesh(end))), x);
dopingFunc = @(x) arrayfun(@(x) dopingFunction(mod(x,mesh(end))),x);
% E0
tempFunc = @(x) arrayfun(@(x) integral(@(y) -coeff * (electronConcentrationFunc(y) - dopingFunc(y)), 0, x),  x);

% tempFunc = @(x) arrayfun(@(x) quadgk(@(y) -coeff*(electronConcentrationFunc(y)-dopingFunc(y)), 0, x),  x);
tic
electricField0 = integral(@(x) tempFunc(x), 0, 1);
toc
disp(['electricField0 cal time:',num2str(toc)]);
% E^h
global VOLTAGE_DROP
electricField = quadgk(@(x) -coeff * (electronConcentration(x)-dopingFunc(x)),x);
electricField = electricField + electricField0 - VOLTAGE_DROP;
electricField = matlabFunction(electricField);

end


%%
function getF(mesh, auxq, electronConcentration, E, n)
global MOBILITY THETA RELAXATION_PARAMETER
F = zeros(N*(n+1),1);
for j = 1:N
    for l = 1:n+1
        [Pl, DPl] = legendreBaseFunction(l-1, mesh(j),mesh(j+1));
        T1 = quadgk(@(x) MOBILITY*E(x) .* electronConcentration(x) .* DPl(x), mesh(j),mesh(j+1));
        T2 = quadgk(@(x) sqrt(THETA * RELAXATION_PARAMETER) * auxq(x) .* DPl(x), mesh(j),mesh(j+1));
        T3 = (MOBILITY * upwindFlux4En(E,electronConcentration, mesh(j+1)) + sqrt(THETA * RELAXATION_PARAMETER) * auxq(mesh(j+1)) ) * Pl(mesh(j+1));
        T4 = (MOBILITY * upwindFlux4En(E,electronConcentration, mesh(j)) + sqrt(THETA * RELAXATION_PARAMETER) * auxq(mesh(j)) ) * Pl(mesh(j));
        F((j-1)*(n+1) + l) = -T1-T2+T3-T4;
    end
end
end
function y = upwindFlux4En(E,electronConcentration,x)
y = max(E(x),0) .* electronConcentration(x) + min(E(x),0) .* electronConcentration(x);
end
