
addpath('scripts');
addpath('functions');
DD = Model('mobilityEqConst.json');
n = 3;
for i = 1:length(N)
    for j = i
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, DD);
        mesh.t = t(j);
        [tempM{i,j}, TV{i,j,5}, e{i,j,5}] = mesh.DDModelDGFunction(2);
    end
end
M{5} = tempM;
clear tempM;