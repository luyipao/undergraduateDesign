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
        initialPriFunction
        degree (1,1) double {mustBeInteger}
        Cells (:,1) Cell
        YR (1,:) double {mustBeNumeric}
        YL (1,:) double {mustBeNumeric}
        CFL = 1;
        epsilon = 1;
        priLegendreFunctions
    end
    methods
        %% generate function
        function obj = Mesh(a, b, N, f, degree)
            if nargin == 0
                ...
            elseif nargin == 5
            f0 = @(x) x;
            f1 = @(x) 1/2 * x.^2;
            f2 = @(x) 1/2 * (x.^3 - x);
            f3 = @(x) 5/8 * x.^4 - 3/4 * x.^2;
            obj.priLegendreFunctions = {f0 f1 f2 f3};
            
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
            %get primitive doping function
            [~,dopingCoeffs] = piecewiseL2Projection(@(x) dopingFunction(x),obj.degree,obj.xa,obj.xb,obj.CellsNum);
			fluxs1 = zeros(N,1);
			fluxs2 = zeros(N,1);
			for j = 1:N
				func  = @(x) 0*x;
				for i = 1:degree+1
					func = @(x) func(x) + dopingCoeffs((j-1)*(degree+1)+i) * sqrt((2*i-1)*obj.meshSize)/2 * obj.priLegendreFunctions{i}((2*x - 2*obj.Xc(j))/obj.meshSize);
				end
				fluxs1(j) = func(obj.X(j+1)) - func(obj.X(j));
				fluxs2(j) = func(obj.X(j));
			end
			fluxs1 = [0;fluxs1];
			fluxs1 = cumsum(fluxs1);
			fluxs1 = fluxs1(1:end-1) - fluxs2;
            dopingCoeffs = reshape(dopingCoeffs,degree+1,length(obj.X)-1);
            obj.initialPriFunction = LegendrePoly(obj.X,dopingCoeffs,degree,fluxs1);
        end
        %% save
        function [] =  save(obj,filename,varargin)
            coeff = obj.coeffs;
            save(filename,'coeff',varargin{:});
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
    end
    methods
        %% deal with DD model
        function obj = DDModelDGFunction(obj)
            Time = [];
            TV = [];
            FLast = obj.getFunc;
            [obj, E, Func] = obj.GK3;
            FNow = Func;
            while sqrt(gaussLegendre(@(x) (FLast(x) - FNow(x)).^2, obj.xa, obj.xb)) > obj.epsilon
                tStart = cputime;
                [obj, E, Func] = obj.GK3;
                tEnd = cputime - tStart;
                Time = [Time tEnd];
                FLast = FNow;
                FNow = Func;
                Y = FNow(obj.Xc);
                TV = [TV sum(abs(Y(2:end)-Y(1:end-1)))];
%                 obj.save('workspace');
%                 save('workspace', 'FNow', 'Time', 'TV', '-append');
            end
            disp('over');
        end
        %% GK
        function [obj, E, Func] = GK3(obj)
            t = obj.CFL * obj.meshSize^2;
            k0 = obj.coeffs;
            x = linspace(0,0.6,1000);
            F = L(obj);
            k1 = k0 + t * F;
            obj.coeffs = k1;
            F = L(obj);
            k2 = 3/4 * k0 + 1/4 * k1 + 1/4 * t * F;
            obj.coeffs = k2;
            [F,E] = L(obj);
            obj.coeffs = 1/3 * k0 + 2/3 * k2 + 2/3 * t * F;
            Func = obj.getFunc;
        end
        % input: obj
        % output: F, electric field E, electron concentration func;
        
        function [F,E] = L(obj)
            n = obj.degree;
            N = obj.CellsNum;
            h = obj.meshSize;
            parforX = obj.X;
            parforXc = obj.Xc;
            parforCoeffs = reshape(obj.coeffs,n+1,N);
            
            CellFuncs = cell(N,1);
            CellValues = zeros(N,4);%[ll lr rl rr]
            
            parforFunc2s = cell(N);
            parforFluxs = zeros(N,1);
            parforCells = obj.Cells;
            parfor (j = 1:N,0)
                % get electron concentration and cell values
                coeffsj = parforCoeffs(:,j);
                Cellj = parforCells(j);
                a = Cellj.a;
                b = Cellj.b;
                fProj = @(x) 0 * x;
                for i = 1:n+1
                    fProj = @(x) fProj(x) + coeffsj(i) * Cellj.basisFunctions{i,1}(x);
                end
                CellFuncs{j} = fProj;
                
                fa = fProj(a);
                fb = fProj(b);
                CellValues(j,:) = [0 fa fb 0];
                % get primitive electron concentration
                % input: obj.coeffs parforX
                % ouput: primitive electron concentration
                
                func2 = @(x) 0*x;
                for i = 1:n+1
                    func2 = @(x) func2(x) + coeffsj(i) * sqrt((2*i-1)*h)/2 * obj.priLegendreFunctions{i}((2*x - 2*parforXc(j))/h);
                end
                parforFluxs(j) = func2(b) - func2(a);
                parforFunc2s{j} = @(x) func2(x);
            end
            % get cell values
            %
            CellValues(2:N-1, 1) = CellValues(1:N-2, 3);
            CellValues(2:N-1, 4) = CellValues(3:N, 2);
            
            % boundary
            CellValues(1, 1) = CellValues(N, 3);
            CellValues(1, 4) = CellValues(2, 2);
            CellValues(N, 1) = CellValues(N-1, 3);
            CellValues(N, 4) = CellValues(1, 2);
            
            
            % get primitive electron concentration
            parforFluxs = [0; parforFluxs];
            parforFluxs = cumsum(parforFluxs);
            for j = 1:N
                parforFluxs(j) = -parforFunc2s{j}(parforX(j)) + parforFluxs(j);
            end
			electronConcentration = LegendrePoly(parforX, parforCoeffs,n,parforFluxs);
			
%             priElectronConcentration = @(x) 0*x;
%             for j = 1:N-1
%                 priElectronConcentration = @(x) priElectronConcentration(x) + (x>=parforX(j) & x < parforX(j+1)) .* (parforFunc2s{j}(x) + parforFluxs(j));
%             end
%             priElectronConcentration = @(x) priElectronConcentration(x) + (x>=parforX(N) & x <= parforX(N+1)) .* (parforFunc2s{N}(x) + parforFluxs(N));
%             
            
            % get electric field
            T1 = gaussLegendre(@(x) electronConcentration.priSolve(x), 0, 0.6)  + gaussLegendre(@(x) electronConcentration.priSolve(x)+electronConcentration.priSolve(0.6),0,0.4);
            T2 = gaussLegendre(@(x) obj.initialPriFunction.priSolve(x), 0, 0.6)  + gaussLegendre(@(x) obj.initialPriFunction.priSolve(x)+obj.initialPriFunction.priSolve(0.6),0,0.4); % set in generate function
            electricField0 = 0.0015 * (T1 - T2);
            E = @(x) -0.0015 * (electronConcentration.priSolve(x) - obj.initialPriFunction.priSolve(x) - electronConcentration.priSolve(0) + obj.initialPriFunction.priSolve(0)) + electricField0 -  1.5;
            
            
            % get auxiliary function
            % input: Cellj.func
            auxCellValues = zeros(N,4); % [ll lr rl rr]
            auxCellFuncs = cell(N,1);
            auxCoeffs = zeros(n+1,N);
            parfor (j = 1:N,0)
                Cellj = parforCells(j);
                temp = zeros(n+1,1);
                a = Cellj.a;
                b = Cellj.b;
                CellValuej = CellValues(j,:);
                for i = 1:n+1
                    tempf = @(x) Cellj.basisFunctions{i,2}(x) .* CellFuncs{j}(x);
                    temp(i) = - gaussLegendre(tempf,a,b);
                end
                
                Fj = CellValuej(4) * Cellj.basisFunctionsBoundaryValues(:,2) - CellValuej(2) * Cellj.basisFunctionsBoundaryValues(:,1) + temp;
                Fj =  0.139219332249189 * Fj;
                auxCoeffs(:,j) = Fj;
                
                auxq = @(x) 0*x;
                for i = 1:n+1
                    auxq = @(x) auxq(x) + Fj(i) * Cellj.basisFunctions{i,1}(x);
                end
                auxCellFuncs{j} = auxq;
                
                auxCellValues(j,:) = [0 auxq(a) auxq(b) 0];
            end
            
            %
            auxCellValues(2:N-1, 1) = auxCellValues(1:N-2, 3);
            auxCellValues(2:N-1, 4) = auxCellValues(3:N, 2);
            % boundary
            auxCellValues(1, 1) = auxCellValues(N, 3);
            auxCellValues(1, 4) = auxCellValues(2, 2);
            auxCellValues(N, 1) = auxCellValues(N-1, 3);
            auxCellValues(N, 4) = auxCellValues(1, 2);
            % input: parforCells,CellValues, auxCellValues
            % L output F very slow
            % temp var: Cellj CellValuej auxCellValuej Fj
            F = zeros(n+1,N);
            for j = 1:N
                Cellj = parforCells(j);
                a = Cellj.a;
                b = Cellj.b;
                CellValuej = CellValues(j,:);
                auxCellValuej = auxCellValues(j,:);
                Eb = E(b);
                T3 = 0.75 * (max(Eb,0) * CellValuej(4) + min(Eb,0) * CellValuej(3));
                T3 = T3 + 0.139219332249189 * auxCellValuej(3);
                Ea = E(a);
                T4 = 0.75 * (max(Ea,0) * CellValuej(2) + min(Ea,0) * CellValuej(1));
                T4 = T4 + 0.139219332249189 * auxCellValuej(1);
                Fj = zeros(n+1,1);
                for l = 1:n+1
                    T1 = gaussLegendre(@(x) 0.75 * E(x) .* CellFuncs{j}(x) .* Cellj.basisFunctions{l,2}(x), a, b);
                    T2 = gaussLegendre(@(x) 0.139219332249189 * auxCellFuncs{j}(x) .* Cellj.basisFunctions{l,2}(x), a, b);
                    Fj(l) = -T1-T2+T3 * Cellj.basisFunctionsBoundaryValues(l,2) - T4 * Cellj.basisFunctionsBoundaryValues(l,1);
                end
                F(:,j) = Fj;
                % free
            end
            F = reshape(F,[],1);
        end
    end
    methods (Access = private)
        % input obj.coeffs
        % output obj.func
        function [Func, obj] = getFunc(obj)
            N = obj.CellsNum;
            n = obj.degree;
            parforCoeffs = reshape(obj.coeffs,n+1,N);
            CellFuncs = cell(N,1);
            parforX = obj.X;
            parforCells = obj.Cells;
            parfor (j = 1:N,0)
                % get electron concentration and cell values
                coeffsj = parforCoeffs(:,j);
                Cellj = parforCells(j);
                fProj = @(x) 0 * x;
                for i = 1:n+1
                    fProj = @(x) fProj(x) + coeffsj(i) * Cellj.basisFunctions{i,1}(x);
                end
                CellFuncs{j} = fProj;
            end
            Func = @(x) 0*x;
            for j = 1:N-1
                Func = @(x) Func(x) + (x>= parforX(j) & x < parforX(j+1)) .* CellFuncs{j}(x);
            end
            Func = @(x) Func(x) + (x>= parforX(N) & x <= parforX(N+1)) .* CellFuncs{N}(x);
            
            obj.func = Func;
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