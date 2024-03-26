function [fProj,fProjCoeffVec]= L2Projection(f,n,a,b)
fProj = @(x) 0 * (x) ;
fProjCoeffVec = zeros(n+1,1);
for i = 0:n
    fBase = legendreBaseFunction(i,a,b);
    fProjCoeffVec(i+1) = quadgk(@(x) fBase(x).*f(x),a,b,'MaxIntervalCount',1e5);
    fProj = @(x) fProj(x) + fProjCoeffVec(i+1) * fBase(x);
end
end

