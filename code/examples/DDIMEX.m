addpath('..\src');
addpath('..\src\functions');

setParameters
a = 0;
b = 0.6;
n = [0 1 2];
timeOrder = 1;
nSteps = zeros(2,5);
N = [25 50 100 200 400 800];
% t = [1.2e-3, 1.8e-3, 2.4e-3, 3.0e-3, 3.6e-3]/3;
t = 1.2e-3;
% M = cell(length(N),length(t));

%% get exact solution
exactMesh = Mesh2(a,b,3200,@(x) dopingFunction(x),2);
exactMesh.t = t;
exactMesh.epsilon = 0.001;
[exactMesh, exactnSteps, exactCOEFFS, exactQCOEFFS] = exactMesh.IMEXGK(3);
exactSolution = exactMesh.getBasisPolys(exactMesh.coeffs);
%%
for k = 1:length(n)
    for i = 1:length(N)
        for j = 1:length(t)
            mesh = Mesh2(a,b,N(i),@(x) dopingFunction(x),n(k));
            mesh.t = t(j);
            mesh.epsilon = 0.001;
            [mesh, nSteps(i,j), COEFFS{timeOrder}{k,i,j}, QCOEFFS{timeOrder}{k,i,j}] = mesh.IMEXGK(timeOrder);
            IMEX{timeOrder}{k,i,j} = mesh;

            % Generate a filename with N and t
            %         electronConcentration = mesh.getBasisPolys(mesh.coeffs);
            %         plot(x,electronConcentration.solve(x));
            %         filename1 = sprintf('DDIMEXRK%dDegree%dmeshCells%dElectronConcentration_t%g.pdf',Order,n, N(i), t(j));
            %         filename1 = fullfile('..\docs\images', filename1);
            %         exportgraphics(gcf, filename1, 'ContentType', 'vector');
            %
            %
            %         % Generate another filename with N and t
            %         E = mesh.getBasisPolys(mesh.Ecoeffs);
            %         plot(x,E.solve(x));
            %         filename2 = sprintf('DDIMEXRK%dDegree%dmeshCells%d_E_t%g.pdf',order, n, N(i), t(j));
            %         filename2 = fullfile('..\docs\images', filename2);
            %         exportgraphics(gcf, filename2, 'ContentType', 'vector');
        end
    end
end

E = zeros(length(N),length(t));
for k = 1:length(n)
    for i = 1:length(N)
        for j = 1:length(t)
            f = IMEX{timeOrder}{k,i,j}.getBasisPolys(IMEX{timeOrder}{k,i,j}.coeffs);
            x = linspace(a,b,1000);
            %         E(i,j) = gaussLegendre(@(x) (f.solve(x)-exactSolution.solve(x)).^2,a , b, 10);
            E(i,j) = sqrt(sum((f.solve(x)-exactSolution.solve(x)).^2))  /  sqrt((2^(i-1)));
        end
    end
    order = (E(1:end-1)-E(2:end))./E(2:end);
    order = log2(order);
    order = num2str(order,'%.4f');
    order = ['------'; order];
    T = table(N', num2str(E, '%.2e'),order );
    filename = sprintf('..\\docs\\tables\\tempMTimeDDIMEXRK%dDegree%d.tex', timeOrder,n(k));
    table2latex(T,filename)
end