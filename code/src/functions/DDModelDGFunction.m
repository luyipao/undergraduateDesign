%output n^h coeff matrix: C; n^h Value matrix electronConcentrationValue;
% space mesh: mesh; Time mesh: T.
function [C, electronConcentrationValue, mesh, T] = DDModelDGFunction(electronConcentrationCoeffVec,auxq, mesh, n, epsilon, CFL)
N = length(mesh) - 1;
meshSize = max(abs(mesh(2:end)-mesh(1:end-1)));
% A
A = zeros((n+1)*N,(n+1)*N);
for j = 1:N
    A_j = zeros(n+1, n+1);
    for k = 1:n+1
        [Pk, DPk] = legendreBaseFunction(k-1,mesh(j),mesh(j+1));
        for l = 1:n+1
            [Pl, DPl] = legendreBaseFunction(l-1,mesh(j),mesh(j+1));
            A_j(k,l) = quadgk(@(x) Pk(x).*Pl(x), mesh(j), mesh(j+1));
        end
    end
    A((j-1)*(n+1)+1 : j*(n+1), (j-1)*(n+1)+1 : j*(n+1)) = A_j;
end

% t
t = CFL * meshSize;

C = [electronConcentrationCoeffVec GK3(electronConcentrationCoeffVec,A,t,mesh,n)] ;
%draw 
[electronConcentration,~,~] = getElectronConcentration(C(:,end), mesh, n);
plot(linspace(0.05,0.6,1000), electronConcentration(linspace(0.05,0.6,1000)));
% T
T = [0];
nLast = getElectronConcentration(C(:,end-1),mesh,n);
nNow = getElectronConcentration(C(:,end),mesh,n);
X = linspace(mesh(1),mesh(end),1000);
electronConcentrationValue = [nLast(X)' nNow(X)'];

while(norm(electronConcentrationValue(:,end-1) - electronConcentrationValue(:,end)) > epsilon)
    T = [T T(end)+t];
    C = [C GK3(C(:,end),A,t,mesh,n)] ;
    nNow = getElectronConcentration(C(:,end),mesh,n);
    electronConcentrationValue = [electronConcentrationValue nNow(X)'];
end

end

function C = GK3(C,A,t,mesh,n)
% global F auxq E electronConcentration priElectronConcentration 
N = length(mesh) - 1;
[electronConcentration,priElectronConcentration,electronConcentrationCells] = getElectronConcentration(C, mesh, n); %ok
auxq = auxiliaryDDModelDGFunction(mesh, electronConcentration, n);%ok
E = getElectricField(priElectronConcentration);% seem ok
F = getF(mesh, auxq, electronConcentration,electronConcentrationCells, E, n);% seem ok
global MOBILITY
for j = 1:N
    B_j = zeros(n+1, n+1);
    for k = 1:n+1
        [Pk, ~] = legendreBaseFunction(k-1,mesh(j),mesh(j+1));
        for l = 1:n+1
            [~, DPl] = legendreBaseFunction(l-1,mesh(j),mesh(j+1));
            B_j(k,l) = quadgk(@(x) MOBILITY * E(x) .* Pk(x).*DPl(x), mesh(j), mesh(j+1));
        end
    end
    B((j-1)*(n+1)+1 : j*(n+1), (j-1)*(n+1)+1 : j*(n+1)) = B_j;
end

k1 = C + t * (A \ F); %- t * A \ (B*C); %%error line k1 can't get correct result.  

[electronConcentration,priElectronConcentration,electronConcentrationCells] = getElectronConcentration(k1, mesh, n);
auxq = auxiliaryDDModelDGFunction(mesh, electronConcentration, n);
E = getElectricField(priElectronConcentration);
F = getF(mesh, auxq, electronConcentration, electronConcentrationCells,E, n);

k2 = 3/4 * C + 1/4 * k1 + 1/4 * t * A \ F;%% k2 already fail.

[electronConcentration,priElectronConcentration,electronConcentrationCells] = getElectronConcentration(k2, mesh, n);
auxq = auxiliaryDDModelDGFunction(mesh, electronConcentration, n);
E = getElectricField(priElectronConcentration);
F = getF(mesh, auxq, electronConcentration, electronConcentrationCells, E, n); %  error 10e8

C = 1/3 * C + 2/3 * k2 + 2/3 * t * A \ F;
end

%%
function electricField = getElectricField(priElectronConcentration)
global ELECTRON_CHARGE DIELECTRIC_PERMITTIVITY
coeff = ELECTRON_CHARGE / DIELECTRIC_PERMITTIVITY;

% E0
T1 = quadgk(@(x) priElectronConcentration(x), 0, 1);
T2 = quadgk(@(x) priDopingFunction(x), 0, 1);
electricField0 = coeff * (T1 - T2);
% E^h
global VOLTAGE_DROP

electricField = @(x) -coeff * (priElectronConcentration(x) - priDopingFunction(x) - priElectronConcentration(0) + priDopingFunction(0)) + electricField0 - VOLTAGE_DROP;
end


%%
function F = getF(mesh, auxq, electronConcentration,electronConcentrationCells, E, n)
global MOBILITY THETA RELAXATION_PARAMETER
N = length(mesh) - 1;
F = zeros(N*(n+1),1);
for j = 1:N
    for l = 1:n+1
        [Pl, DPl] = legendreBaseFunction(l-1, mesh(j),mesh(j+1));
        T1 = quadgk(@(x) MOBILITY * E(x) .* electronConcentration(x) .* DPl(x), mesh(j),mesh(j+1));
        T2 = quadgk(@(x) sqrt(THETA * RELAXATION_PARAMETER) * auxq(x) .* DPl(x), mesh(j),mesh(j+1));
        
        if j == N
            upwindFlux = max(E(mesh(j+1)), 0) .* electronConcentrationCells{1}(mesh(1)) + min(E(mesh(j+1)), 0) .* electronConcentrationCells{j}(mesh(j+1));
        else
            upwindFlux = max(E(mesh(j+1)),0) .* electronConcentrationCells{j+1}(mesh(j+1)) + min(E(mesh(j+1)), 0) .* electronConcentrationCells{j}(mesh(j+1));
        end
        T3 = (MOBILITY * upwindFlux  + sqrt(THETA * RELAXATION_PARAMETER) * auxq(mesh(j+1))) * Pl(mesh(j+1));
        
        if j == 1
            upwindFlux = max(E(mesh(j)),0) .* electronConcentrationCells{j}(mesh(j)) + min(E(mesh(j)), 0) .* electronConcentrationCells{end}(mesh(end));
            T4 = (MOBILITY * upwindFlux + sqrt(THETA * RELAXATION_PARAMETER) * auxq(mesh(j))) * Pl(mesh(j));
        else
            upwindFlux = max(E(mesh(j)),0) .* electronConcentrationCells{j}(mesh(j)) + min(E(mesh(j)), 0) .* electronConcentrationCells{j-1}(mesh(j));
            T4 = (MOBILITY * upwindFlux  + sqrt(THETA * RELAXATION_PARAMETER) * auxq(mesh(j))) * Pl(mesh(j));
        end
        F((j-1)*(n+1) + l) = -T1-T2+T3-T4;
    end
end
end

