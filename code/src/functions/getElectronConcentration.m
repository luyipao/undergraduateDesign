%% input: electronConcentrationCoeff C; mesh; poly max degree n
% function type: [) 
function [electronConcentration,priElectronConcentration,electronConcentrationCells] = getElectronConcentration(C,mesh,n)
N = length(mesh) - 1;
electronConcentration = @(x) 0*x;
priElectronConcentration = @(x) 0*x;
electronConcentrationCells = cell(N,1);
flux = 0;
func2 = @(x) 0*x;

for j = 1:N
    func1 = @(x) 0*x;
    for i = 1:n+1
        [Pi,~] = legendreBaseFunction(i-1, mesh(j),mesh(j+1));
        func1 = @(x) func1(x) + C((j-1)*(n+1)+i) * Pi(x);
    end
    func2 = matlabFunction(int(sym(func1)));
    flux = flux - func2(mesh(j));
    if j==1
        flux = 0;
    end
    electronConcentrationCells{j} = func1;
    if j==N
        electronConcentration = @(x) electronConcentration(x) + (x>=mesh(j) & x <= mesh(j+1)) .* func1(x);
        priElectronConcentration = @(x) priElectronConcentration(x) + (x>=mesh(j) & x <= mesh(j+1)) .* (func2(x) + flux);
    else
        electronConcentration = @(x) electronConcentration(x) + (x>=mesh(j) & x < mesh(j+1)) .* func1(x);
        priElectronConcentration = @(x) priElectronConcentration(x) + (x>=mesh(j) & x < mesh(j+1)) .* (func2(x) + flux);
    end
    flux = flux + func2(mesh(j+1));
end
end