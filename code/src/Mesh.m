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
        func
        auxFunction
        degree (1,1) double {mustBeInteger}
        Cells (:,1) Cell
        YR (1,:) double {mustBeNumeric}
        YL (1,:) double {mustBeNumeric}
        CFL = 0.1429;
        epsilon = 0.1;
    end
    methods
        %% generate function
        function obj = Mesh(a, b, N, f, degree)
            if nargin == 0
                ...
            elseif nargin == 5
            obj.xa = a;
            obj.xb = b;
            obj.meshSize = (obj.xb-obj.xa) / N;
            obj.CellsNum = N;
            obj.degree = degree;
            obj.X = linspace(a, b, N+1);
            obj.Xc = (obj.X(1:end-1) + obj.X(2:end)) / 2;
            obj.coeffs = [];
%             obj.func = @(x) 0 * x;
            for i = 1:N
                obj.Cells(i) = Cell(obj.X(i),obj.X(i+1),degree);
                [~,fProjCoeffVec] = L2Projection(f,degree,obj.X(i),obj.X(i+1));
                obj.Cells(i).coeffs = fProjCoeffVec;
                obj.coeffs = [obj.coeffs; fProjCoeffVec];
%                 if i==N
%                     obj.func = @(x) obj.func(x) + (x>=obj.X(i) & x<=obj.X(i+1)) .* fProj(x);
%                 else
%                     obj.func = @(x) obj.func(x) + (x>=obj.X(i) & x<obj.X(i+1)) .* fProj(x);
%                 end
            end
            else
                error('Invalid number of inputs.')
            end
            obj = obj.getFunc;
            obj = obj.getNodesValues;
            obj = obj.slopeLimiter;
            obj = obj.getFunc;
            obj = obj.getNodesValues;
            % obj = obj.ENOreconstruction;
            %obj = obj.auxiliaryDDModelDGFunction;
        end
        %% ENO reconstruction
        % accuracy: 2m+1
        % input: mesh
        % output: reconstructed function value at mesh nodes in both sides
        function obj = ENOreconstruction(obj)
            m = round(obj.degree / 2) + 1;
            for j = 1:obj.CellsNum
                pPos = obj.auxENOSolver(obj.func,j,m);
                pNeg = obj.auxENOSolver(obj.func,j+1,m);
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
            dividedDiffX = obj.Xc(2:end) - obj.Xc(1:end-1);
            if i < 1
                L = cumsum(dividedDiffX(end+i:end),'reverse');
                if j > N
                    R = cumsum(dividedDiffX(1:j-N));
                    mesh = [obj.Xc(1)-L obj.Xc obj.Xc(end)+R];
                else
                    mesh = [obj.Xc(1)-L obj.Xc(1:j)];
                end
            else
                if j > N
                    R = cumsum(dividedDiffX(1:j-N));
                    mesh = [obj.Xc(i:end) obj.Xc(end)+R];
                else
                    mesh = obj.Xc(i:j);
                end
            end
        end
        %% set Coeffs
    end
    methods (Access = public)
        %% getFunc
        % input: coeff
        % output: func in Cell and mesh
        function obj = getFunc(obj)
            obj.func = @(x) 0*x;
            for j = 1:obj.CellsNum
                fProj = @(x) 0 * x;
                Cellj = obj.Cells(j);
                for i = 1:obj.degree+1
                    fProj = @(x) fProj(x) + obj.coeffs((j-1)*(obj.degree+1) + i) * Cellj.basisFunctions{i,1}(x);
                end
                obj.Cells(j).func = @(x) fProj(x);
                if j == obj.CellsNum
                    obj.func = @(x) obj.func(x) + (x >= Cellj.a & x <= Cellj.b) .* fProj(x);
                else
                    obj.func = @(x) obj.func(x) + (x >= Cellj.a & x < Cellj.b) .* fProj(x);
                end
            end
        end
        %% getNodesValues
        function obj = getNodesValues(obj)
            for j = 1:obj.CellsNum
                fProj = @(x) obj.Cells(j).func(x);
                obj.Cells(j).rl = fProj(obj.X(j+1));
                obj.Cells(j).lr = fProj(obj.X(j));
                if j == 1
                    obj.Cells(obj.CellsNum).rr = obj.Cells(j).lr;
                    obj.Cells(j+1).ll = obj.Cells(j).rl;
                elseif j == obj.CellsNum
                    obj.Cells(1).ll = obj.Cells(j).rl;
                    obj.Cells(j-1).rr = obj.Cells(j).lr;
                else
                    obj.Cells(j+1).ll = obj.Cells(j).rl;
                    obj.Cells(j-1).rr = obj.Cells(j).lr;
                end

            end
        end
        %% auxENOsolver
        function p = auxENOSolver(obj,f,j,m)
            % init
            f = @(x) f(mod(x-obj.xa, obj.xb-obj.xa) + obj.xa);  % extend f periodly to all x bias
            a = zeros(2*m+1,1);
            b = zeros(2*m+1,1);
            c = zeros(2*m+1,1);
            Q = cell(2*m+1,1);
            kmin = j;
            kmax = j;
            if j > obj.CellsNum
                Q{1} = @(x) f(obj.Xc(j-obj.CellsNum)) + 0*x;
            else
                Q{1} = @(x) f(obj.Xc(j)) + 0*x;
            end
            % inductively
            for i = 2:2*m+1
                % nth divided differences of f
                ax = obj.getLocalMesh(kmin,kmax+1);
                ay = f(ax);
                a(i) = obj.dividedDiff(ax,ay);
                bx = obj.getLocalMesh(kmin-1,kmax);
                by = f(bx);
                b(i) = obj.dividedDiff(bx,by);
                if abs(a(i))+1 >= abs(b(i))  % error  when a(i) should equalt to b(i), small noisy may cause b(i) > a(i)
                    c(i) = b(i);
                    kmin = kmin-1;
                else
                    c(i) = a(i);
                    kmax = kmax+1;
                end
                % product
                product = @(x) 1 + 0 * x;
                for xk = obj.getLocalMesh(kmin,kmax)
                    product = @(x) product(x) .* (x - xk);
                end
                % Q
                Q{i} = @(x) Q{i-1}(x) + c(i) *  product(x);
            end
            % output
            p = @(x) Q{end}(x);
%             plot(linspace(obj.X(j),obj.X(j+1)),p(linspace(obj.X(j),obj.X(j+1))),linspace(obj.X(j),obj.X(j+1)),obj.func(linspace(obj.X(j),obj.X(j+1))),'--');
%             legend('reconstructed','origin');
        end
        %% dividedDiff
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
        %% getElectronConcentration
        % input:
        % output:
        function priElectronConcentration = getElectronConcentration(obj)
            % init
            N = obj.CellsNum;
            n = obj.degree;
            % output
            priElectronConcentration = @(x) 0*x;
            flux = 0;
            for j = 1:N
                syms x
                func1 = @(x) 0*x;
                Cellj = obj.Cells(j);
                for i = 1:n+1
                    func1 = @(x) func1(x) + obj.coeffs((j-1)*(n+1)+i) * Cellj.basisFunctions{i,1}(x);
                end
                func2 = matlabFunction(int(sym(func1)),'vars', x);
                flux = flux - func2(obj.X(j));
                if j==1
                    flux = 0;
                end
                if j==N
                    priElectronConcentration = @(x) priElectronConcentration(x) + (x>=obj.X(j) & x <= obj.X(j+1)) .* (func2(x) + flux);
                else
                    priElectronConcentration = @(x) priElectronConcentration(x) + (x>=obj.X(j) & x < obj.X(j+1)) .* (func2(x) + flux);
                end
                flux = flux + func2(obj.X(j+1));
            end
        end
        %% ElectricField
        function electricField = getElectricField(priElectronConcentration)
            coeff = 0.0015;
            % E0
            T1 = gaussLegendre(@(x) priElectronConcentration(x), 0, 1);
            T2 = gaussLegendre(@(x) priDopingFunction(x), 0, 1);
            electricField0 = coeff * (T1 - T2);
            % E^h
            electricField = @(x) -coeff * (priElectronConcentration(x) - priDopingFunction(x) - priElectronConcentration(0) + priDopingFunction(0)) + electricField0 -  1.5;
        end
    end
    methods
        %% deal with DD model
        function obj = DDModelDGFunction(obj)
            
            obj = obj.GK3;
        end
        %% GK
        function obj = GK3(obj)
            t = obj.CFL * obj.meshSize;
            k0 = obj.coeffs;
            priElectronConcentration = obj.getElectronConcentration();
            E = getElectricField(priElectronConcentration);% seem ok
            F = getF(obj,E);
            k1 = k0 + t * F;
            obj.coeffs = k1;
            obj = obj.getFunc;
            obj = obj.getNodesValues;
            %obj = obj.ENOreconstruction;
            obj = obj.auxiliaryDDModelDGFunction;
            priElectronConcentration = obj.getElectronConcentration();
            E = getElectricField(priElectronConcentration);% seem ok
            F = getF(obj,E);
            obj = obj.getFunc;
            obj = obj.getNodesValues;
            k2 = 3/4 * k0 + 1/4 * k1 + 1/4 * t * F;
            obj.coeffs = k2;
            obj = obj.getFunc;
            obj = obj.getNodesValues;
            %obj = obj.ENOreconstruction;
            obj = obj.auxiliaryDDModelDGFunction;
            priElectronConcentration = obj.getElectronConcentration();
            E = getElectricField(priElectronConcentration);% seem ok
            F = getF(obj,E);
            obj.coeffs = 1/3 * k0 + 2/3 * k2 + 2/3 * t * F;
            obj = obj.getFunc;
            %obj = obj.ENOreconstruction;
            obj = obj.auxiliaryDDModelDGFunction;
        end
        %% internal function getF
        function F = getF(obj,E)
            N = obj.CellsNum;
            n = obj.degree;
            mesh = obj.X;
            F = zeros(N*(n+1),1);
            for j = 1:N
                Cellj = obj.Cells(j);
                for l = 1:n+1
                    T1 = gaussLegendre(@(x) 0.75 * E(x) .* Cellj.func(x) .* Cellj.basisFunctions{l,2}(x), mesh(j),mesh(j+1));
                    T2 = gaussLegendre(@(x) 0.139219332249189 * Cellj.auxFunction(x) .* Cellj.basisFunctions{l,2}(x), mesh(j),mesh(j+1));
                    
                    T3 = 0.75 * (max(E(mesh(j+1)),0) * Cellj.rr + min(E(mesh(j+1)),0) * Cellj.rl);
                    T3 = T3 + 0.139219332249189 * Cellj.auxrl;
                    T3 = T3 * Cellj.basisFunctionsBoundaryValues(l,2);
                    
                    T4 = 0.75 * (max(E(mesh(j+1)),0)*Cellj.lr + min(E(mesh(j)),0)*Cellj.ll);
                    T4 = T4 + 0.139219332249189 * Cellj.auxll;
                    T4 = T4 * Cellj.basisFunctionsBoundaryValues(l,1);
                    
                    F((j-1)*(n+1) + l) = -T1-T2+T3-T4;
                end
            end
        end
        function obj = auxiliaryDDModelDGFunction(obj)
            n = obj.degree;
            F = zeros(obj.CellsNum*(n+1),1);
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                % Fj = Cellj.rr * Cellj.basisFunctionsBoundaryValues(:,2) - Cellj.lr * Cellj.basisFunctionsBoundaryValues(:,1);
                Fj = obj.func(obj.X(j+1)) * Cellj.basisFunctionsBoundaryValues(:,2) - obj.func(obj.X(j)) * Cellj.basisFunctionsBoundaryValues(:,1);
                for i = 1:n+1
                    tempf = @(x) Cellj.basisFunctions{i,2}(x) .* Cellj.func(x);
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
            YY = zeros(1,1000);
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                Y = Y + (XX>=Cellj.a & XX<Cellj.b) .* Cellj.auxFunction(XX);
                YY = YY + (XX>=Cellj.a & XX<Cellj.b) .* Cellj.func(XX);
            end
            Y = Y + (XX == obj.xb) .* Cellj.auxFunction(XX);
            YY  = YY + (XX == obj.xb) .* Cellj.func(XX);
            plot(XX,Y,XX,YY);
            legend('auxiliray','aimed function');
        end
        % limiter
        % change coeffs not function
        function obj = slopeLimiter(obj)
            localAverages = zeros(obj.CellsNum,1);
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                localAverages(j) = gaussLegendre(Cellj.func, Cellj.a, Cellj.b) / obj.meshSize;
            end
            M = 6.553417322323334e+08 * obj.meshSize^2;
            N = obj.CellsNum;
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                ur = Cellj.rl - localAverages(j);
                ul = localAverages(j) - Cellj.lr;
                urmod = localAverages(j) + obj.minmod([ul localAverages(mod(j,N)+1)-localAverages(j) localAverages(j)-localAverages(mod(j-2+N,N)+1)], M);
                ulmod =  localAverages(j) - obj.minmod([ur localAverages(mod(j,N)+1)-localAverages(j) localAverages(j)-localAverages(mod(j-2+N,N)+1)], M);
                if abs(urmod - Cellj.rl) < 1 && abs(ulmod - Cellj.lr) < 1
                    continue;
                end
                Cj = zeros(obj.degree+1,1);
                Cj(1) = sqrt(obj.meshSize) * localAverages(j);
                Cj(2) = (urmod - ulmod) / (2*sqrt( 3/obj.meshSize ));
                Cj(3) = urmod - Cj(1) / sqrt(obj.meshSize) - Cj(2) * sqrt(3/obj.meshSize);
                obj.Cells(j).coeffs = Cj;
                obj.coeffs((j-1)*(obj.degree+1)+1:j*(obj.degree+1)) = Cj; 
            end
        end
        function y = minmod(~,v,M)
        if abs(v(1)) <= M
            y = v(1);
        else
            y = all(sign(v) == sign(v(1))) * sign(v(1)) * min(abs(v));
        end
end
    end
end