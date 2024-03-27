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
A = PiPjIntegral;
I = eye(N);
A = kron(I,A);
% cal F
F = zeros(N*(n+1),1);
for j = 1:N
    for l = 1:n+1
        [~,DPl] = legendreBaseFunction(l-1,mesh(j),mesh(j+1));
        index = (j-1)*(n+1) + l;
        F(index) = nFunc(mesh(j+1)) * basisFunctionValue(l,2) - nFunc(mesh(j)) * basisFunctionValue(l,2) - quadgk(@(x) nFunc(x) .* DPl(x), mesh(j), mesh(j+1));
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




