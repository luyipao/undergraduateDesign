% N: cells number
function [fPiecewiseProj,fPiecewiseProjCoeffVec] = piecewiseL2Projection(f,n,a,b,N)
X = linspace(a,b,N+1);
fPiecewiseProj = @(x) 0 * x;
fPiecewiseProjCoeffVec = zeros(N*(n+1),1);
for i = 1:N
    [fProj, fProjCoeffVec] = L2Projection(f,n,X(i),X(i+1));
    fPiecewiseProj = @(x) fPiecewiseProj(x) + (x>=X(i) & x<X(i+1)) .* fProj(x);
    fPiecewiseProjCoeffVec((i-1)*(n+1)+1:i*(n+1)) = fProjCoeffVec;
end
fPiecewiseProj = @(x) fPiecewiseProj(x) + (x == b) * f(b);
end
