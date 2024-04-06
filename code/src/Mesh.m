classdef Mesh
    properties
        xa (1,1) double {mustBeNumeric, mustBeReal}
        xb (1,1) double {mustBeNumeric, mustBeReal}
        meshSize (1,1) double
        CellsNum (1,1) double {mustBeInteger, mustBeFinite}
        X (1,:) double {mustBeNumeric}
        Xc (1,:) double
        coeffs (:,1) double {mustBeNumeric}
        auxCoeffs (:,1) double
        auxFunction (:,1) double
        degree (1,1) double {mustBeInteger}
        Cells (:,1) Cell
        YR (1,:) double {mustBeNumeric}
        YL (1,:) double {mustBeNumeric}
        CFL = 0.1429;
        epsilon = 0.1;
    end
    methods
        %% generate function
        function obj = Mesh(a, b, N, degree)
            if nargin == 0
                ...
            elseif nargin == 4
            obj.xa = a;
            obj.xb = b;
            obj.meshSize = (obj.xb-obj.xa) / N;
            obj.CellsNum = N;
            obj.degree = degree;
            obj.X = linspace(a, b, N+1);
            obj.Xc = (obj.X(1:end-1) + obj.X(2:end)) / 2;
            obj.coeffs = zeros(N*(degree+1),1);
            for i = 1:N
                obj.Cells(i) = Cell(obj.X(i),obj.X(i+1),degree);
            end
            else
                error('Invalid number of inputs.')
            end
        end
        %% ENO reconstruction
        % accuracy: 2m+1
        % input: mesh; aimed function f
        % output: reconstructed function value at mesh nodes in both sides
        function obj = ENOreconstruction(obj,f)
            m = round(obj.degree / 2);
            for j = 1:obj.CellsNum
                pPos = obj.auxENOSolver(f,j,m);
                pNeg = obj.auxENOSolver(f,j+1,m);
                x = obj.X(j+1);
                obj.Cells(j).rr = pPos(x);
                obj.Cells(j).rl = pNeg(x);
                if j == obj.CellsNum
                    obj.Cells(1).lr = pPos(x);
                    obj.Cells(1).ll = pNeg(x);
                else
                    obj.Cells(j+1).lr = pPos(x);
                    obj.Cells(j+1).ll = pNeg(x);
                end
            end
        end
        
        %% get Cell Info
        function Cell = getCell(obj, i)
            idx = mod(i - 1 + obj.CellsNum, obj.CellsNum) + 1;
            Cell = obj.Cells(idx);
        end
        function x = getNode(obj,i)
            i = mod(i - 1 + obj.CellsNum, obj.CellsNum) + 1;
            x = obj.X(i);
        end
        
        %% get extended mesh nodes
        function mesh = getLocalMesh(obj,i,j)
            N = obj.CellsNum;
            dividedDiffX = obj.X(2:end) - obj.X(1:end-1);
            if i < 1
                L = cumsum(dividedDiffX(end+i:end),'reverse');
                if j > N + 1
                    R = cumsum(dividedDiffX(1:j-N - 1));
                    mesh = [obj.X(1)-L obj.X obj.X(end)+R];
                else
                    mesh = [obj.X(1)-L obj.X(1:j)];
                end
            else
                if j > N + 1
                    R = cumsum(dividedDiffX(1:j-N - 1));
                    mesh = [obj.X(i:end) obj.X(end)+R];
                else
                    mesh = obj.X(i:j);
                end
            end
        end
        %% set Coeffs
        function coeffs = setCoeffs(obj, coeffs)
            ...
        end
    end
    methods (Access = private)
        %% auxENOsolver
        function p = auxENOSolver(obj,f,j,m)
            % init
            a = zeros(2*m+1,1);
            b = zeros(2*m+1,1);
            c = zeros(2*m+1,1);
            Q = cell(2*m+1,1);
            kmin = j;
            kmax = j;
            Q{1} = @(x) f(x);
            % inductively
            for i = 2:2*m+1
                % nth divided differences of f
                ax = obj.getLocalMesh(kmin,kmax+1);
                ay = f(ax);
                a(i) = obj.dividedDiff(ax,ay);
                bx = obj.getLocalMesh(kmin-1,kmax+1);
                by = f(bx);
                b(i) = obj.dividedDiff(bx,by);
                if abs(a(i)) >= abs(b(i))
                    c(i) = b(i);
                    kmin = kmin-1;
                else
                    c(i) = a(i);
                    kmax = kmax+1;
                end
                % product
                product = @(x) 0 * x;
                for xk = obj.getLocalMesh(kmin,kmax)
                    product = @(x) product(x) .* (x - xk);
                end
                % Q
                Q{i} = @(x) Q{i-1}(x) + c(i) *  product(x);
            end
            % output
            p = @(x) Q{end}(x);
        end
        %%
        function result = dividedDiff(~, x, y)
            n = length(x);
            diff_table = zeros(n, n);
            diff_table(:, 1) = y(:);
            
            for j = 2:n
                for i = j:n
                    diff_table(i, j) = (diff_table(i, j - 1) - diff_table(i - 1, j - 1)) / (x(i) - x(i - j + 1));
                end
            end
            result = diff_table(n, n);
        end
        %%
        % input: mesh, aimed function f
        % output: auxiliary function q no coeffs vector both in mesh and
        % cell
    end
    methods
        %% deal with DD model
        function obj = DDModelDGFunction(obj,f)
            t = obj.CFL * obj.meshSize;
            [~, obj.coeffs] = piecewiseL2Projection(f,obj.X,obj.degree);
            [electronConcentration,priElectronConcentration,electronConcentrationCells] = getElectronConcentration(obj.coeffs, obj.X, obj.degree);
            E = getElectricField(priElectronConcentration);% seem ok
        end
        function F = getF(obj,E,electronConcentration,electronConcentrationCells)
            N = obj.CellsNum;
            n = obj.degree;
            mesh = obj.X;
            for j = 1:N
                for l = 1:n+1
                    [Pl, DPl] = legendreBaseFunction(l-1, mesh(j),mesh(j+1));
                    T1 = gaussLegendre(@(x) 0.75 * E(x) .* electronConcentration(x) .* DPl(x), mesh(j),mesh(j+1));
                    T2 = gaussLegendre(@(x) 0.139219332249189 * auxq(x) .* DPl(x), mesh(j),mesh(j+1));
                    Cellj = obj.Cells(j);
                    if j == N
                        upwindFlux = max(E(mesh(j+1)), 0) .* electronConcentrationCells{1}(mesh(1)) + min(E(mesh(j+1)), 0) .* electronConcentrationCells{j}(mesh(j+1));
                    else
                        upwindFlux = max(E(mesh(j+1)),0) .* electronConcentrationCells{j+1}(mesh(j+1)) + min(E(mesh(j+1)), 0) .* electronConcentrationCells{j}(mesh(j+1));
                    end
                    T3 = (MOBILITY * upwindFlux  + sqrt(THETA * RELAXATION_PARAMETER) * auxq(mesh(j+1))) * Pl(mesh(j+1));
                    
                    if j == 1
                        upwindFlux = max(E(mesh(j)),0) .* electronConcentrationCells{j}(mesh(j)) + min(E(mesh(j)), 0) .* electronConcentrationCells{end}(mesh(end));
                        T4 = (MOBILITY * upwindFlux + sqrt(THETA * RELAXATION_PARAMETER) * auxq(mesh(j))) * Pl(mesh(j));
                    else
                        upwindFlux = max(E(mesh(j)),0) .* electronConcentrationCells{j}(mesh(j)) + min(E(mesh(j)), 0) .* electronConcentrationCells{j-1}(mesh(j));
                        T4 = (MOBILITY * upwindFlux  + sqrt(THETA * RELAXATION_PARAMETER) * auxq(mesh(j))) * Pl(mesh(j));
                    end
                    F((j-1)*(n+1) + l) = -T1-T2+T3-T4;
                end
            end
            
        end
        
        function obj = auxiliaryDDModelDGFunction(obj,f)
            n = obj.degree;
            F = zeros(obj.CellsNum*(n+1),1);
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                Fj = Cellj.rr * Cellj.basisFunctionsBoundaryValues(:,2) - Cellj.lr * Cellj.basisFunctionsBoundaryValues(:,1);
                for i = 1:n+1
                    tempf = @(x) Cellj.basisFunctions{i,2}(x) .* f(x);
                    Fj(i) = Fj(i) - gaussLegendre(tempf,Cellj.a,Cellj.b);
                end
                Fj =  0.139219332249189 * Fj;
                Cellj.auxCoeffs = Fj;
                F((j-1)*(n+1)+1 : j*(n+1)) = Fj;
                % cal Cell value
                auxq = @(x) 0*x;
                for i = 1:numel(Cellj.basisFunctions(:,1))
                    auxq = @(x) auxq(x) + Fj(i) * Cellj.basisFunctions{i,1}(x);
                end
                Cellj.auxFunction =@(x) auxq(x);
                Cellj.auxrl = auxq(obj.X(j+1));
                Cellj.auxlr = auxq(obj.X(j));
                if j==obj.CellsNum
                    obj.Cells(1).auxll = Cellj.auxrl;
                    obj.Cells(j-1).auxrr = Cellj.auxlr;
                elseif j==1
                    obj.Cells(obj.CellsNum).auxrr = Cellj.auxlr;
                    obj.Cells(j+1).auxll = Cellj.auxrl;
                else
                    obj.Cells(j+1).auxll = Cellj.auxrl;
                    obj.Cells(j-1).auxrr = Cellj.auxlr;
                end
                obj.Cells(j) = Cellj;
            end
            % get Mesh values
            obj.auxCoeffs = F;
        end
        function draw(obj)
            XX = linspace(obj.xa,obj.xb,1000);
            Y = zeros(1,1000);
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                Y = Y + (XX>Cellj.a & XX<Cellj.b) .* Cellj.auxFunction(XX);
            end
            plot(XX,Y);
            legend('auxiliray');
        end
    end
end