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
        initialCoeffs (:,1) double
        func
        auxFunction
        initialFunction
        degree (1,1) double {mustBeInteger}
        Cells (:,1) Cell
        YR (1,:) double {mustBeNumeric}
        YL (1,:) double {mustBeNumeric}
        CFL = 1;
        epsilon = 100;
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
            obj.CFL = 1/(2*degree+1);
            obj.coeffs = [];
            for i = 1:N
                obj.Cells(i) = Cell(obj.X(i),obj.X(i+1),degree);
                [fProj,fProjCoeffVec] = L2Projection(f,degree,obj.X(i),obj.X(i+1));
                obj.Cells(i).coeffs = fProjCoeffVec;
                obj.coeffs = [obj.coeffs; fProjCoeffVec];
            end
            else
                error('Invalid number of inputs.')
            end
            obj = obj.getFunc;
            obj = obj.getNodesValues;
            obj = obj.auxiliaryDDModelDGFunction;
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
    end
    methods (Access = public)
        %% getFunc
        % input: coeff
        % output: func in Cell and mesh
        function obj = getFunc(obj)
            obj.func = @(x) 0*x;
            ParforCells = obj.Cells;
            ParforCoeffs = obj.coeffs;
            n = obj.degree;
            N = obj.CellsNum;
            parfor j = 1:obj.CellsNum
                fProj = @(x) 0 * x;
                for i = 1:n+1
                    fProj = @(x) fProj(x) + ParforCoeffs((j-1)*(n+1) + i) * ParforCells(j).basisFunctions{i,1}(x);
                end
                ParforCells(j).func = @(x) fProj(x);
            end
            obj.Cells = ParforCells;
            parforFunc = @(x) 0*x;
            for j = 1:N-1
                    parforFunc = @(x) parforFunc(x) + (x >= obj.X(j) & x < obj.X(j+1)) .* ParforCells(j).func(x);
            end
            parforFunc = @(x) parforFunc(x) + (x >= obj.X(N) & x <= obj.X(N+1)) .* ParforCells(N).func(x);
            obj.func = @(x) parforFunc(x);
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
            f0 = @(x) x;
            f1 = @(x) 1/2 * x.^2;
            f2 = @(x) 1/2 * (x.^3 - x);
            f3 = @(x) 5/8 * x.^4 - 3/4 * x.^2;
            priLegendreFunctions = {f0 f1 f2 f3};
            priElectronConcentration = @(x) 0*x;
            
            parforCoeffs = obj.coeffs;
            h = obj.meshSize;
            N = obj.CellsNum;
            n = obj.degree; 
            parforXc = obj.Xc;
            parforX = obj.X;
            
            parforFunc2s = cells(N);
            parforFluxs = zeros(N,1);
            parfor j = 1:N
                func2 = @(x) 0*x;
                for i = 1:n+1
                    func2 = @(x) func2(x) +parforCoeffs((j-1)*(n+1)+i) * sqrt((2*i-1)*h)/2 * priLegendreFunctions{i}((2*x - 2*parforXc(j))/h);
                end
                parforFluxs(j) = func2(parforX(j+1)) - func2(parforX(j));
                parforFunc2s{j} = @(x) func2(x);
            end
            parforFluxs = [0; parforFluxs];
            parforFluxs = cumsum(parforFluxs);
            parfor j = 1:N
                parforFluxs(j) = -parforFunc2s{j}(parforX(j)) + parforFluxs(j);
            end
            for j = 1:N-1
                priElectronConcentration = @(x) priElectronConcentration(x) + (x>=parforX(j) & x < parforX(j+1)) .* (parforFunc2s{j}(x) + parforFluxs(j));
            end
            priElectronConcentration = @(x) priElectronConcentration(x) + (x>=parforX(N) & x <= parforX(N+1)) .* (parforFunc2s{N}(x) + parforFluxs(N));

        end
    end
    methods
        %% deal with DD model
        function obj = DDModelDGFunction(obj)
            C(:,1) = obj.coeffs;
            F{1} = obj.func;
            obj = obj.GK3;
            C(:,2) = obj.coeffs;
            F{2} = obj.func;
            while sqrt(gaussLegendre(@(x) (F{end-1}(x) - F{end}(x)).^2, obj.xa, obj.xb)) > obj.epsilon
                tStart = cputime;
                obj = obj.GK3;
                tEnd = cputime - tStart;
                disp(tEnd);
                C(:,end+1) = obj.coeffs;
                F{end+1} = obj.func;
            end
            disp('over');
        end
        %% GK
        function obj = GK3(obj)
            t = obj.CFL * obj.meshSize^2;
            k0 = obj.coeffs;
            x = linspace(0,0.6,1000);
            priElectronConcentration = obj.getElectronConcentration;% run slowlt
            E = obj.getElectricField(priElectronConcentration);% seem ok
            F = getF(obj,E);
            k1 = k0 + t * F;
            obj.coeffs = k1;
            obj = obj.getFunc;
            obj = obj.getNodesValues;
            obj = obj.auxiliaryDDModelDGFunction;
            priElectronConcentration = obj.getElectronConcentration;
            E = getElectricField(priElectronConcentration);% seem ok
            F = getF(obj,E);
            obj = obj.getFunc;
            obj = obj.getNodesValues;
            k2 = 3/4 * k0 + 1/4 * k1 + 1/4 * t * F;
            obj.coeffs = k2;
            obj = obj.getFunc;
            obj = obj.getNodesValues;
            obj = obj.auxiliaryDDModelDGFunction;
            priElectronConcentration = obj.getElectronConcentration();
            E = getElectricField(priElectronConcentration);% seem ok
            F = getF(obj,E);
            obj.coeffs = 1/3 * k0 + 2/3 * k2 + 2/3 * t * F;
            obj = obj.getFunc;
            obj = obj.getNodesValues;
            obj = obj.auxiliaryDDModelDGFunction;
        end
        %% ElectricField pro
        function electricField = getElectricField(obj,priElectronConcentration)
            f0 = @(x) x;
            f1 = @(x) 1/2 * x.^2;
            f2 = @(x) 1/2 * (x.^3 - x);
            f3 = @(x) 5/8 * x.^4 - 3/4 * x.^2;
            priLegendreFunctions = {f0 f1 f2 f3};
            % E0
            n=obj.degree;
            N=obj.CellsNum;
            T1 = gaussLegendre(@(x) priElectronConcentration(x), 0, 1);
            [dopingProj,dopingcCoeffs] = piecewiseL2Projection(@(x) dopingFunction(x),obj.degree,obj.xa,obj.xb,obj.CellsNum);
            priDopingFunction = @(x) 0*x;
            flux = 0;
            for j = 1:N
                func2 = @(x) 0*x;
                for i = 1:n+1
                    func2 = @(x) func2(x) + dopingcCoeffs((j-1)*(n+1)+i) * sqrt((2*i-1)*obj.meshSize)/2 * priLegendreFunctions{i}((2*x - 2*obj.Xc(j))/obj.meshSize);
                end
                flux = flux - func2(obj.X(j));
                if j==1
                    flux = 0;
                end
                if j==N
                    priDopingFunction = @(x) priDopingFunction(x) + (x>=obj.X(j) & x <= obj.X(j+1)) .* (func2(x) + flux);
                else
                    priDopingFunction = @(x) priDopingFunction(x) + (x>=obj.X(j) & x < obj.X(j+1)) .* (func2(x) + flux);
                end
                flux = flux + func2(obj.X(j+1));
            end
            T2 = gaussLegendre(@(x) priDopingFunction(x), 0, 1);
            electricField0 = 0.0015 * (T1 - T2);
            electricField = @(x) -0.0015 * (priElectronConcentration(x) - priDopingFunction(x) - priElectronConcentration(0) + priDopingFunction(0)) + electricField0 -  1.5;
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
        %% auxiliaryDDModelDGFunctionx
        % output: auxfunction Nodes values; auxfunction
        function obj = auxiliaryDDModelDGFunction(obj)
            n = obj.degree;
            F = zeros(obj.CellsNum*(n+1),1);
            % get coeff
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                Fj = Cellj.rr * Cellj.basisFunctionsBoundaryValues(:,2) - Cellj.lr * Cellj.basisFunctionsBoundaryValues(:,1);
                for i = 1:n+1
                    tempf = @(x) Cellj.basisFunctions{i,2}(x) .* Cellj.func(x);
                    Fj(i) = Fj(i) - gaussLegendre(tempf,Cellj.a,Cellj.b);
                end
                Fj =  0.139219332249189 * Fj;
                Cellj.auxCoeffs = Fj;
                F((j-1)*(n+1)+1 : j*(n+1)) = Fj;
                obj.Cells(j) = Cellj;
            end
            %cal Cell value
            obj = obj.setAuxFunction;
            % get Mesh values
            obj.auxCoeffs = F;
        end
        %input: aux coeffs
        % output: auxfunction and its nodes value;
        function obj = setAuxFunction(obj)
            obj.auxFunction = @(x) 0*x;
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                auxq = @(x) 0*x;
                for i = 1:numel(Cellj.basisFunctions(:,1))
                    auxq = @(x) auxq(x) + Cellj.auxCoeffs(i) * Cellj.basisFunctions{i,1}(x);
                end
                Cellj.auxFunction =@(x) auxq(x);
                obj.auxFunction = @(x) obj.auxFunction(x) + (Cellj.a <= x & x < Cellj.b) .* auxq(x);
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
        end
        % draw the image of function and auxFunction
    end
    methods (Access = private)
        %% ENO reconstruction
        % accuracy: 2m+1
        % input: mesh
        % output: reconstructed function value at mesh nodes in both sides
        function obj = ENOreconstruction(obj)
            m = round(obj.degree / 2) + 1;
            L = obj.xb - obj.xa;
            %            f = @(x) obj.func(mod(x-obj.xa,L)+obj.xa);
            f = @(x) obj.func(mod(x,L));
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
        function obj = auxENOreconstruction(obj)
            m = round(obj.degree/2)+1;
            L = obj.xb - obj.xa;
            %            f = @(x) obj.auxFunction(mod(x-obj.xa,L)+obj.xa);
            f = @(x) obj.auxFunction(mod(x,L));
            for j = 1:obj.CellsNum
                pPos = obj.auxENOSolver(f,j,m);
                pNeg = obj.auxENOSolver(f,j+1,m);
                x = obj.X(j+1);
                obj.Cells(j).auxrr = pPos(x);
                obj.Cells(j).auxrl = pNeg(x);
                if j == obj.CellsNum
                    obj.Cells(1).auxlr = pPos(x);
                    obj.Cells(1).auxll = pNeg(x);
                else
                    obj.Cells(j+1).auxlr = pPos(x);
                    obj.Cells(j+1).auxll = pNeg(x);
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
            kmin = zeros(2*m+1,1);
            kmax = zeros(2*m+1,1);
            kmin(1) = j;
            kmax(1) = j;
            if j == obj.CellsNum + 1
                Q{1} = @(x) f(obj.Xc(1)) + 0*x;
            else
                Q{1} = @(x) f(obj.Xc(j)) + 0*x;
            end
            % inductively
            for i = 2:2*m+1
                % nth divided differences of f
                ax = obj.getLocalMesh(kmin(i-1),kmax(i-1)+1);
                ay = f(ax);
                a(i) = obj.dividedDiff(ax,ay);
                bx = obj.getLocalMesh(kmin(i-1)-1,kmax(i-1));
                by = f(bx);
                b(i) = obj.dividedDiff(bx,by);
                if abs(a(i))+1 >= abs(b(i))  % error  when a(i) should equalt to b(i), small noisy may cause b(i) > a(i)
                    c(i) = b(i);
                    kmin(i) = kmin(i-1)-1;
                    kmax(i) = kmax(i-1);
                else
                    c(i) = a(i);
                    kmin(i) = kmin(i-1);
                    kmax(i) = kmax(i-1)+1;
                end
                % product
                product = @(x) 1 + 0 * x;
                for xk = obj.getLocalMesh(kmin(i-1),kmax(i-1))
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
        
        %% get extended mesh nodes
        % outputcondition: mesh may exceed [xa,xb]
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
                urmod = localAverages(j) + obj.minmod([ur localAverages(mod(j,N)+1)-localAverages(j) localAverages(j)-localAverages(mod(j-2+N,N)+1)], M);
                ulmod =  localAverages(j) - obj.minmod([ul localAverages(mod(j,N)+1)-localAverages(j) localAverages(j)-localAverages(mod(j-2+N,N)+1)], M);
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
        % input: auxfunction
        % output: auxcoeff
        function obj = auxSlopeLimiter(obj)
            localAverages = zeros(obj.CellsNum,1);
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                localAverages(j) = gaussLegendre(Cellj.auxFunction, Cellj.a, Cellj.b) / obj.meshSize;
            end
            M = 5.731745892056191e+09 * obj.meshSize^2;
            N = obj.CellsNum;
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                ur = Cellj.auxrl - localAverages(j);
                ul = localAverages(j) - Cellj.auxlr;
                urmod = localAverages(j) + obj.minmod([ur localAverages(mod(j,N)+1)-localAverages(j) localAverages(j)-localAverages(mod(j-2+N,N)+1)], M);
                ulmod =  localAverages(j) - obj.minmod([ul localAverages(mod(j,N)+1)-localAverages(j) localAverages(j)-localAverages(mod(j-2+N,N)+1)], M);
                if abs(urmod - Cellj.auxrl) < 1 && abs(ulmod - Cellj.auxlr) < 1
                    continue;
                end
                Cj = zeros(obj.degree+1,1);
                Cj(1) = sqrt(obj.meshSize) * localAverages(j);
                Cj(2) = (urmod - ulmod) / (2*sqrt( 3/obj.meshSize ));
                Cj(3) = urmod - Cj(1) / sqrt(obj.meshSize) - Cj(2) * sqrt(3/obj.meshSize);
                obj.Cells(j).auxCoeffs = Cj;
                obj.auxCoeffs((j-1)*(obj.degree+1)+1:j*(obj.degree+1)) = Cj;
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
    
    methods
        function draw(obj)
            XX = linspace(obj.xa,obj.xb,1000);
            Y = zeros(1,1000);
            YY = zeros(1,1000);
            for j = 1:obj.CellsNum
                Cellj = obj.Cells(j);
                Y = Y + (XX>=Cellj.a & XX<Cellj.b) .* Cellj.auxFunction(XX);
                YY = YY + (XX>=Cellj.a & XX<Cellj.b) .* Cellj.func(XX);
            end
            Y = Y + (XX == obj.xb) .* obj.Cells(end).auxFunction(XX);
            YY  = YY + (XX == obj.xb) .* obj.Cells(end).func(XX);
            hold on;
            plot(XX,Y,'--');
            plot(XX,YY);
            hold off;
            legend('auxiliray','aimed function');
        end
        function drawError(obj)
            plot(linspace(0,0.6,1000),dopingFunction(linspace(0,0.6,1000)) - obj.func(linspace(0,0.6,1000)));
        end
    end
end