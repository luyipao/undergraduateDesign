E = zeros(length(N),length(t));
for i = 1:length(N)
    for j = 1:length(t)
        f = IMEX2M{i,j}.getBasisPolys(IMEX2M{i,j}.coeffs);
        x = linspace(a,b,1000);
        E(i,j) = sqrt(sum((f.solve(x)-exactSolution.solve(x)).^2)./sum(exactSolution.solve(x).^2));
    end
end
order = (E(1:end-1))./E(2:end);
order = log2(order);
order = num2str(order,'%.2f');
order = ['----'; order] ;
T = table(N', num2str(E, '%.2e'),order );
filename = sprintf('..\\docs\\tables\\DDIMEXRK2Degree%d.tex', n);
table2latex(T,filename)

%% draw time when stable reached
% T1 = table(num2str(t', '%.1e'), nSteps(1,:)',num2str(t' .* nSteps(1,:)', '%.4f'), 'VariableNames', {'$\Delta t$', 'nt', 't'});
% filename1 = sprintf('..\\docs\\tables\\DDIMEXRK3Degree%dN100T.tex', n);
% T2 = table(num2str(t', '%.1e'), nSteps(2,:)',num2str(t' .* nSteps(2,:)', '%.4f'), 'VariableNames', {'$\Delta t$', 'nt', 't'});
% filename2 = sprintf('..\\docs\\tables\\DDIMEXRK3Degree%dN200T.tex', n);
% table2latex(T1,filename1);
% table2latex(T2,filename2);

%% TVDRK VS IMEX RK
% x = linspace(0, 0.6);
% % N = 100
% f1 = IMEXM{4,1}.getBasisPolys(IMEXM{4,1}.coeffs);
% E1 = IMEXM{4,1}.getBasisPolys(IMEXM{4,1}.Ecoeffs);
% 
% % N = 200
% 
% 
% f2 = IMEXM{3,1}.getBasisPolys(IMEXM{3,1}.coeffs);
% E2 = IMEXM{3,1}.getBasisPolys(IMEXM{3,1}.Ecoeffs);
% 
% plot(x,f1.solve(x), x, f2.solve(x), 'o', 'linewidth', 1);
% lgd = legend('N=100', 'N=200', 'Location', 'north');
% yticks(0:100000:600000);
% setFormat(lgd, 'x', 'n');
% filename = sprintf('IMEXNn.pdf');
% filename = fullfile('..\..\论文\docs\figure',filename);
% exportgraphics(gcf, filename, 'ContentType', 'vector');
% 
% plot(x,E1.solve(x), x, E2.solve(x), 'o', 'linewidth', 1);
% lgd = legend('N=100', 'N=200', 'Location', 'north');
% 
% setFormat(lgd, 'x', 'E');
% filename = sprintf('IMEXNE.pdf');
% filename = fullfile('..\..\论文\docs\figure',filename);
% exportgraphics(gcf, filename, 'ContentType', 'vector');
% 
% function setFormat(lgd, xLabelName, yLabelName)
% lgd.FontSize = 11; % 
% lgd.Location = 'north'; %
% set(gca,'XMinorTick', 'on', 'YMinorTick', 'on', 'linewidth',1.5); % change axes width
% 
% set(gca,'yticklabel',get(gca,'ytick'));
% set(gca,'xticklabel',get(gca,'xtick'));
% xlabel(xLabelName);
% ylabel(yLabelName);
% end