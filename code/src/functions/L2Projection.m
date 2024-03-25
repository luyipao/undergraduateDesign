function fProj = L2Projection(f,n,a,b)
fProj = @(x) 0 * (x) ;
for i = 0:n
    fi = legendreBaseFunction(i,a,b);
    fProj = @(x) fProj(x) + quadgk(@(x) fi(x).*f(x),a,b) * fi(x);
end
end

