
addpath('C:\Users\12192\OneDrive - zju.edu.cn\documents\homeworks\undergraduateDesign\code\src\functions');
setParameters
a = 0;
b = 0.6;
n = 2;
N = 100;

mesh = Mesh(a,b,N,@(x) dopingFunction(x),n);
mesh = mesh.DDModelDGFunction();


