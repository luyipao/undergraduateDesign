%% getElectronConcentration
% input:
% output:
function priElectronConcentration = getElectronConcentration(obj)
% init
N = obj.CellsNum;
n = obj.degree;
% output
priElectronConcentration = @(x) 0*x;
flux = 0;
for j = 1:N
    syms x
    func1 = @(x) 0*x;
    Cellj = obj.Cells(j);
    for i = 1:n+1
        func1 = @(x) func1(x) + obj.coeffs((j-1)*(n+1)+i) * Cellj.basisFunctions{i,1}(x);
    end
    func2 = matlabFunction(int(sym(func1)),'vars', x);
    flux = flux - func2(obj.X(j));
    if j==1
        flux = 0;
    end
    if j==N
        priElectronConcentration = @(x) priElectronConcentration(x) + (x>=obj.X(j) & x <= obj.X(j+1)) .* (func2(x) + flux);
    else
        priElectronConcentration = @(x) priElectronConcentration(x) + (x>=obj.X(j) & x < obj.X(j+1)) .* (func2(x) + flux);
    end
    flux = flux + func2(obj.X(j+1));
end
end
function priElectronConcentration = getElectronConcentration(obj)
f0 = @(x) x;
f1 = @(x) 1/2 * x.^2;
f2 = @(x) 1/2 * (x.^3 - x);
f3 = @(x) 5/8 * x.^4 - 3/4 * x.^2;
priLegendreFunctions = {f0 f1 f2 f3};
N = obj.CellsNum;
n = obj.degree;
priElectronConcentration = @(x) 0*x;
flux = 0;
for j = 1:N
    func2 = @(x) 0*x;
    for i = 1:n+1
        func2 = @(x) func2(x) + obj.coeffs((j-1)*(n+1)+i) * sqrt((2*i-1)*obj.meshSize)/2 * priLegendreFunctions{i}((2*x - 2*obj.Xc(j))/obj.meshSize);
    end
    flux = flux - func2(obj.X(j));
    if j==1
        flux = 0;
    end
    if j==N
        priElectronConcentration = @(x) priElectronConcentration(x) + (x>=obj.X(j) & x <= obj.X(j+1)) .* (func2(x) + flux);
    else
        priElectronConcentration = @(x) priElectronConcentration(x) + (x>=obj.X(j) & x < obj.X(j+1)) .* (func2(x) + flux);
    end
    flux = flux + func2(obj.X(j+1));
end