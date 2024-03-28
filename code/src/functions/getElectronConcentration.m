%% input: electronConcentrationCoeff C; mesh; poly max degree n
function electronConcentration = getElectronConcentration(C,mesh,n)
N = length(mesh) - 1;
electronConcentration = @(x) 0*x;
for j = 1:N
    for i = 1:n+1
        [Pi,~] = legendreBaseFunction(i-1, mesh(j),mesh(j+1));
        electronConcentration = @(x) electronConcentration(x) + C((j-1)*(n+1)+i) * (x>=mesh(j) & x < mesh(j+1)) .* Pi(x);
    end
end
electronConcentration = @(x) electronConcentration(x) + (x == mesh(end)) * electronConcentration(mesh(1));
end