addpath('..\src');
addpath('..\src\functions');

setParameters
a = 0;
b = 0.6;
n = 2;
nSteps = zeros(2,5);
N = [25 50 100 200 400 800];
%t = [1.2e-3, 1.8e-3, 2.4e-3, 3.0e-3, 3.6e-3]/3;
t = 1.2e-3;
M = cell(length(N),length(t));
% get exact solution
Fmesh = Mesh2(a,b,3200,@(x) dopingFunction(x),n);
Fmesh.t = t;
[Fmesh, exactnSteps] = Fmesh.IMEXGK(3);
exactSolution = Fmesh.getBasisPolys(Fmesh.coeffs);
for i = 1:length(N)
    for j = 1:length(t)
        mesh = Mesh2(a,b,N(i),@(x) dopingFunction(x),n);
        mesh.t = t(j);
        [mesh, nSteps(i,j)] = mesh.IMEXGK(3);
        M{i,j} = mesh;
        x = linspace(0,0.6,1000);
        plot(x,electronConcentration.solve(x));
        
        % Generate a filename with N and t
        filename1 = sprintf('DDIMEXRK3Degree3meshCells%dElectronConcentration_t%g.pdf', N(i), t(j));
        print('-dpdf',fullfile('..\docs\images', filename1));
        
        E = mesh.getBasisPolys(mesh.Ecoeffs);
        plot(x,E.solve(x));
        
        % Generate another filename with N and t
        filename2 = sprintf('DDIMEXRK3Degree3meshCells%d_E_t%g.pdf', N(i), t(j));
        print('-dpdf',fullfile('..\docs\images', filename2));
    end
end
% E = zeros(length(N),length(t));
% for i = 1:length(N)
%     for j = 1:length(t)
%         f = M{i,j}.getBasisPolys(M{i,j}.coeffs);
%         x = linspace(a,b,1000);
%         E(i,j) = sqrt(sum((f.solve(x)-exactSolution.solve(x)).^2)/1000);
%     end
% end

E = zeros(length(N),length(t));
for i = 1:length(N)
    for j = 1:length(t)
        f = M{i,j}.getBasisPolys(M{i,j}.coeffs);
        x = linspace(a,b,1000);
        E(i,j) = sqrt(sum((f.solve(x)-exactSolution.solve(x)).^2)./sum(exactSolution.solve(x).^2));
    end
end
T = table(N, num2str(E, '%.2e'), num2str(order,'%.2f'));
filename = sprintf('..\\docs\\tables\\DDIMEXRK3Degree%d.tex', 0);
table2latex(T,filename)