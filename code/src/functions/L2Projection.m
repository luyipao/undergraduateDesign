function [fProj,fProjCoeffVec]= L2Projection(f,n,a,b)
fProj = @(x) 0 * (x) ;
fProjCoeffVec = zeros(n+1,1);
for i = 0:n
    fBasis = legendreBaseFunction(i,a,b);
    fProjCoeffVec(i+1) = gaussLegendre(@(x) fBasis(x).*f(x),a,b);
    fProj = @(x) fProj(x) + fProjCoeffVec(i+1) * fBasis(x);
end
end

