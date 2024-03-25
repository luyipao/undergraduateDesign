% N: [a,b]分的段数

function fPiecewiseProj = piecewiseL2Projection(f,n,a,b,N)
X = linspace(a,b,N+1);
fPiecewiseProj = @(x) 0 * x;
for i = 1:N
    fProj = L2Projection(f,n,X(i),X(i+1));
    fPiecewiseProj = @(x) fPiecewiseProj(x) + (x>=X(i) & x<X(i+1)) .* fProj(x);
end
fPiecewiseProj = @(x) fPiecewiseProj(x) + (x == b) * f(b);
end
