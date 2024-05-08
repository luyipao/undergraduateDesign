
addpath('scripts');
addpath('functions');
setParameters
a = 0;
b = 0.6;
n = 2;
nSteps = zeros(2,5);
N = [100 200];
t = [1.6e-5 4.2e-6];
nmodel = 3;
DD = Model('mobilityEqConst.json');
for i = 1:length(N)
    for j = i
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, DD);
        mesh.t = t(j);
        [tempM{i,j}, TV{i,j,1}, e{i,j,1}] = mesh.DDModelDGFunction(2);
    end
end
M{1} = tempM;
clear tempM;

N = 100;

DD2 = Model('mobilityDpDopingFunction.json');
for i = 1:length(N)
    for j = i
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, DD2);
        mesh.t = t(j);
        [tempM{i,j}, TV{i,j,2}, e{i,j,2}] = mesh.DDModelDGFunction(2);
    end
end
M{2} = tempM;
clear tempM;

HF = Model('HFmobilityEqConst.json');
for i = 1:length(N)
    for j = i
        mesh = Mesh(a,b,N(i),@(x) dopingFunction(x),n, HF);
        mesh.t = t(j);
        [tempM{i,j}, TV{i,j,3}, e{i,j,3}] = mesh.DDModelDGFunction(2);
    end
end
M{3} = tempM;
clear tempM;

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
