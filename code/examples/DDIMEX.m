addpath('..\src');
setParameters
a = 0;
b = 0.6;
n = 3;
nSteps = zeros(2,5);
N = [100 200];
t = [1.2e-3, 1.8e-3, 2.4e-3, 3.0e-3, 3.6e-3]/3;
for i = 1:2
    for j = 1:5
        mesh = Mesh2(a,b,N(i),@(x) dopingFunction(x),n);
        mesh.t = t(j);
        [mesh, nSteps(i,j)] = mesh.IMEXGK(3);
        x = linspace(0,0.6,1000);
        electronConcentration = mesh.getBasisPolys(mesh.coeffs);
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