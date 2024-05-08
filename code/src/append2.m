
addpath('scripts');
addpath('functions');
DD = Model('mobilityEqConst.json');
for i = 1:length(N)
    for j = i
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, DD);
        mesh.t = t(j);
        [tempM{i,j}, TV{i,j,4}, e{i,j,4}] = mesh.DDModelDGFunction(2);
    end
end
M{4} = tempM;
clear tempM;