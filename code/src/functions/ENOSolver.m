
f+u+n+ction [posP, negP] = ENOSolver(posf,negf,alpha,j,m,mesh)

posP = auxENOSolver(posf,alpha,j,m,mesh);
negP = auxENOSolver(negf,alpha,j+1,m,mesh);
x = mesh(j+1);
posFluxf = posP(x);
for k = 1:m-1
    
end

function p = auxENOSolver(f,alpha,j,m,mesh)
%% init
f = @(x) 1/2 * (f(x) + alpha*x);
a = zeros(2*m+1,1);
b = zeros(2*m+1,1);
c = zeros(2*m+1,1);
Q = cell(2*m+1,1);
kmin = j;
kmax = j;
Q{0} = @(x) f(x);
%% inductively
for i = 2:2*m+1
    % nth divided differences of f
    ax = mesh(kmin:kmax+1);
    ay = f(x);
    a(i) = dividedDiff(ax,ay);
    bx = mesh(kmin-1:kmax+1);
    by = f(y);
    b(i) = dividedDiff(bx,by);
    if abs(a(i)) >= abs(b(i))
        c(i) = b(i);
        kmin = kmin-1;
    else
        c(i) = a(i);
        kmax = kmax+1;
    end
    % product
    product = @(x) 0 * x;
    for xk = mesh(kmin,kmax)
        product = @(x) product .* (x - xk);
    end
    Q{i} = @(x) Q{i-1}(x) + c(i) *  product;
end
%% output
p = @(x) Q{end}(x);
end

function result = dividedDiff(x, y)
if length(x) == 1
    result = y;
else
    result = dividedDiff(x(2:end),y(2:end)) - dividedDiff(x(1:end-1),y(1:end-1));
    result = result / (x(end) - x(1));
end
end
