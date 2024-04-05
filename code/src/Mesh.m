classdef Mesh
    properties
        xa (1,1) double {mustBeNumeric, mustBeReal}
        xb (1,1) double {mustBeNumeric, mustBeReal}
        meshSize (1,1) double
        CellsNum (1,1) double {mustBeInteger, mustBeFinite}
        X (1,:) double {mustBeNumeric}
        Xc (1,:) double
        coeffs (:,1) double {mustBeNumeric}
        degree (1,1) double {mustBeInteger}
        Cells (:,1) Cell
        YR (1,:) double {mustBeNumeric}
        YL (1,:) double {mustBeNumeric}
    end
    methods
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
                obj.Cells(i) = Cell();
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
            m = round(obj.n / 2);
            for j = 1:obj.N
                pPos = auxENOSolver(f,j,m);
                pNeg = auxENOSolver(f,j+1,m);
                x = obj.X(j+1);
                obj.Cells(j).rr = pPos(x);
                obj.Cells(j).rl = pNeg(x);
                if j == obj.N
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
            Q{0} = @(x) f(x);
            % inductively
            for i = 2:2*m+1
                % nth divided differences of f
                ax = obj.getLocalMesh(kmin:kmax+1);
                ay = f(x);
                a(i) = dividedDiff(ax,ay);
                bx = obj.getLocalMesh(kmin-1:kmax+1);
                by = f(y);
                b(i) = dividedDiff(bx,by);
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
                    product = @(x) product .* (x - xk);
                end
                % Q
                Q{i} = @(x) Q{i-1}(x) + c(i) *  product;
            end
            % output
            p = @(x) Q{end}(x);
        end
        %%
        function result = dividedDiff(x, y)
            if length(x) == 1
                result = y;
            else
                result = dividedDiff(x(2:end),y(2:end)) - dividedDiff(x(1:end-1),y(1:end-1));
                result = result / (x(end) - x(1));
            end
        end
    end
end