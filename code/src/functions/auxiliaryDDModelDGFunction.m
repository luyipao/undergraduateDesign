% input: mesh; poly max degree n; nFunc
% output: auxiliary fucntion q coeff;


function [auxq,auxiliaryVarCoeffVec] = auxiliaryDDModelDGFunction(mesh, nFunc, n)
%% initlize
global RELAXATION_PARAMETER THETA
N = length(mesh) - 1;
meshSize = max(abs(mesh(2:end)-mesh(1:end-1)));
[basisFunctionValue, PiDPjIntegral, PiPjIntegral] = getLegendreBasisInfo(n,meshSize);

%% AC=F, C is auxiliaryVarCoeffVec
% cal A
% 预先分配矩阵
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
% cal F
F = zeros(N*(n+1),1);
for j = 1:N
    for l = 1:n+1
        [Pl,DPl] = legendreBaseFunction(l-1,mesh(j),mesh(j+1));
        index = (j-1)*(n+1) + l;
        F(index) = nFunc(mesh(j+1)) * Pl(mesh(j+1)) - nFunc(mesh(j)) * Pl(mesh(j)) - quadgk(@(x) nFunc(x) .* DPl(x), mesh(j), mesh(j+1));
        F(index) = F(index) * sqrt(RELAXATION_PARAMETER * THETA);
    end
end
% cal C
auxiliaryVarCoeffVec = A\F;

%% output auxiliary function q
auxq = @(x) 0 * x;

for j = 1:N
    if j == 1
        for l = 1:n+1
            [Pl,~] = legendreBaseFunction(l-1,mesh(j),mesh(j+1));
            auxq = @(x) auxq(x) + auxiliaryVarCoeffVec((j-1)*(n+1)+l) * (x >= mesh(j) & x <= mesh(j+1)) .* Pl(x);
        end
    else
        for l = 1:n+1
            [Pl,~] = legendreBaseFunction(l-1,mesh(j),mesh(j+1));
            auxq = @(x) auxq(x) + auxiliaryVarCoeffVec((j-1)*(n+1)+l) *(x > mesh(j) & x <= mesh(j+1)) .* Pl(x);
        end
    end
end
end




