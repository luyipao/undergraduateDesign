% %% IMEX 1
% % get aux coeffs
timeOrder = 3;
[~, I] = size(exactCOEFFS);
A = exactMesh.meshSize * kron(eye(exactMesh.CellsNum),exactMesh.massMatrix);
A = diag(1./diag(A));
A = sparse(A);
for ii = 1:I
    exactQCOEFFS(:,ii) = A*exactMesh.BPos*exactCOEFFS(:,ii);
end
for k = 1:length(n)
    E = zeros(length(N), length(t));
    for i = 1:length(N)
        for j = 1:length(t)
            tempMesh = IMEX{timeOrder}{k,i,j};
            [~, I] = size(COEFFS{timeOrder}{k,i,j});
            A = tempMesh.meshSize * kron(eye(tempMesh.CellsNum),tempMesh.massMatrix);
            A = diag(1./diag(A));
            A = sparse(A);
            for ii = 1:I
                QCOEFFS{timeOrder}{k,i,j}(:,ii) = A*tempMesh.BPos*COEFFS{timeOrder}{k,i,j}(:,ii);
            end
        end
    end
end
% consider LInf in [0, T] rather than T point
for k = 1:length(n)
    E = zeros(length(N), length(t));
    for i = 1:length(N)
        for j = 1:length(t)
            res = 99999999999999999999;
            SUM = 0;
            for iii = 1:250
                f = IMEX{timeOrder}{k,i,j}.getBasisPolys(COEFFS{timeOrder}{k,i,j}(:,iii));
                exactf = exactMesh.getBasisPolys(exactCOEFFS(:,iii));
                g =  IMEX{timeOrder}{k,i,j}.getBasisPolys(numq{timeOrder}{k,i,j}(:,iii));
                exactg = exactMesh.getBasisPolys(q(:,iii));
                x = linspace(a,b,10000);
                temp = (f.solve(x)-exactf.solve(x)).^2;

                temp = sqrt(sum(temp));
                SUM = SUM + (sum((g.solve(x)- exactg.solve(x)).^2));
                % draw
%                 x = linspace(a, b, 1000);
%                 figure(1)
%                 plot(x, exactg.solve(x), x, g.solve(x) , '--');
%                 figure(2)
%                 plot(x, exactf.solve(x), x, f.solve(x), '--');
                if temp < res
                    res = temp;
                end
            end
            SUM = sqrt(SUM * t);
            E(i,j) = res + SUM;
        end
    end
%         C = (2^k * E(end) - E(end-1)) / (2^k-1);
%     order = (E(1:end-1) - C) ./ (E(2:end) - C);
    order = E(1:end-1) ./ E(2:end);
    order = log(order) / log(2);
    order = num2str(order,'%.4f');
    order = ['------'; order];
    T = table(N', num2str(E, '%.2e'),order );
    filename = sprintf('..\\docs\\tables\\LInfDDIMEXRK%dDegree%d.tex', timeOrder,n(k));
    table2latex(T,filename)
end

%% IMEX 2 & 3
%  n is degree of polynomial; N is mesh num; timeorder is IMEX order;
%  output: order of IMEX cal
% for k = 1:length(n)
%     E = zeros(length(N),length(t));
%     for i = 1:length(N)
%         for j = 1:length(t)
%             [~, I] = size(COEFFS{timeOrder}{k,i,j});
%             x = linspace(a,b,10000);
%             res = 99999999999999999999;
%             constTemp = 0;
%             for ii = 1:I
%                 f = IMEX{timeOrder}{k,i,j}.getBasisPolys(COEFFS{timeOrder}{k,i,j}(:,iii));
%                 exactf = exactMesh.getBasisPolys(exactCOEFFS(:,iii));
%                 temp = sqrt(sum((f.solve(x)-exactf.solve(x)).^2));
%                 if constTemp < sqrt(sum((exactf.solve(x)).^2))
%                     constTemp = sqrt(sum((exactf.solve(x)).^2));
%                 end
%                 if temp < res
%                     res = temp;
%                 end
%             end
%             E(i,j) = res / constTemp;
%         end
%     end
%     
%     %      C = mean((2^k * E(2:end) - E(1:end-1)) ./ (2^k - 1));
%     %     C = (2^k * E(end) - E(end-1)) / (2^k-1);
%     %     order = (E(1:end-1) - C) ./ (E(2:end) - C);
%     order = E(1:end-1) ./ E(2:end);
%     order = log(order) / log(2)
%     order = num2str(order,'%.8f');
%     order = ['----------'; order];
%     T = table(N', num2str(E, '%.2e'),order );
%     filename = sprintf('..\\docs\\tables\\tempMTimeDDIMEXRK%dDegree%d.tex', timeOrder,n(k));
%     table2latex(T,filename)
% end

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