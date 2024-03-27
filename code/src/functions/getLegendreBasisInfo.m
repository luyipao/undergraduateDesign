% input: poly max degree n; interval length h;
% output: basis function value at boundary points. basis function integral.

function [basisFunctionValue, basisFunctionPiDPjIntegral, basisFunctionPiPjIntegral] = getLegendreBasisInfo(n,h)
basisFunctionValue = zeros(n+1,2);
basisFunctionPiDPjIntegral = zeros(n+1,n+1);
basisFunctionPiPjIntegral = zeros(n+1,n+1);
for i = 1:n+1
    [Pi,~] = legendreBaseFunction(i-1,-h/2,h/2);
    basisFunctionValue(i,:) = [Pi(-h/2),Pi(h/2)];
    for j = 1:n+1
        [Pj, DPj] = legendreBaseFunction(j-1,-h/2,h/2);
        basisFunctionPiDPjIntegral(i,j) = quadgk(@(x) Pi(x).*DPj(x), -h/2,h/2);
        basisFunctionPiPjIntegral(i,j) = quadgk(@(x) Pi(x).*Pj(x), -h/2,h/2);
    end
end
end