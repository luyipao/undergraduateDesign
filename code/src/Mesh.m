classdef Mesh
    properties
        xa (1,1) double {mustBeNumeric, mustBeReal}
        xb (1,1) double {mustBeNumeric, mustBeReal}
        degree (1,1) double {mustBeInteger}
        meshSize (1,1) double
        CellsNum (1,1) double {mustBeInteger, mustBeFinite}
        X (1,:) double {mustBeNumeric}
        Xc (1,:) double
        coeffs (:,1) double {mustBeNumeric}
        auxCoeffs (:,1) double
        Ecoeffs	(:,1) double
        initialCoeffs (:,1) double
        func
        auxFunction
        initialFunction
        
        classicalLegendrePolys
        priLegendreFunctions
        diffBasisFuncs
        basisBoundaryValues
        PnDPm
        PPDP
        Cells (:,1) Cell
        CFL = 0.2;
        epsilon = 0.01;
        t
        
        semiModel Model
    end
    methods
        %% generate function
        function obj = Mesh(a, b, N, f, degree, model)
            if nargin == 0
                ...
            elseif nargin >= 5
            f0 = @(x) x;
            f1 = @(x) 1/2 * x.^2;
            f2 = @(x) 1/2 * (x.^3 - x);
            f3 = @(x) 5/8 * x.^4 - 3/4 * x.^2;
            obj.priLegendreFunctions = {f0 f1 f2 f3};
            f0 = @(x) 0*x + 1;
            f1 = @(x) x;
            f2 = @(x) 1.5 * x.^2 - 0.5;
            f3 = @(x) 2.5 * x.^3 - 1.5 * x;
            obj.classicalLegendrePolys = {f0 f1 f2 f3};
            obj.classicalLegendrePolys = obj.classicalLegendrePolys(1:degree+1);
            obj.xa = a;
            obj.xb = b;
            obj.meshSize = (obj.xb-obj.xa) / N;
            obj.t = obj.CFL * obj.meshSize^2;
            temp = 2 / obj.meshSize * [sqrt(1/obj.meshSize) sqrt(3/obj.meshSize) sqrt(5/obj.meshSize) sqrt(7/obj.meshSize)];
            obj.diffBasisFuncs = {@(x) 0*x; @(x) 0*x + temp(2); @(x) 3*temp(3)*x; @(x) temp(4) *(7.5*x.^2 - 1.5)};
            obj.diffBasisFuncs = obj.diffBasisFuncs(1:degree+1);
            
            obj.CellsNum = N;
            obj.degree = degree;
            obj.X = linspace(a, b, N+1);
            obj.Xc = (obj.X(1:end-1) + obj.X(2:end)) / 2;
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
            if nargin == 6
                obj.semiModel = model;
            end
            %get primitive doping function
            [~,dopingCoeffs] = piecewiseL2Projection(@(x) dopingFunction(x),obj.degree,obj.xa,obj.xb,obj.CellsNum);
            obj.initialCoeffs = dopingCoeffs;
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
            
            obj.initialFunction = LegendrePoly(obj.X,dopingCoeffs,degree,fluxs1);
            % get integration of Pn and DPm
            % Pn is scaled standardized classical Legendre polynomials
            obj.PnDPm = zeros(degree+1,degree+1);
            for i = 1:degree+1
                for k = 0:degree+1
                    idx = 1 + k * 2;
                    if (i+idx <= degree+1)
                        obj.PnDPm(i,i+idx) = 2 * sqrt((2*i-1)*(2*i+2*idx-1)) / obj.meshSize;
                    else
                        break;
                    end
                end
            end
            % cell basis function in boundary value;
            % basis functions refer scaled Orthonormal Legendre polynomails
            obj.basisBoundaryValues = [sqrt(1/obj.meshSize) sqrt(3/obj.meshSize) sqrt(5/obj.meshSize) sqrt(7/obj.meshSize)];
            obj.basisBoundaryValues = obj.basisBoundaryValues(1:degree+1);
            obj.basisBoundaryValues = [obj.basisBoundaryValues .* (-1).^(0:degree); obj.basisBoundaryValues]';
            
            % scaled orthonormal legendre polys integral
            syms x
            A = zeros(obj.degree+1, obj.degree+1, obj.degree+1);
            P = legendreP(0:obj.degree+1,x);
            DP = diff(P);
            for i = 1:obj.degree+1
                for j = 1:obj.degree+1
                    for k = 1:obj.degree+1
                        f = P(i) * P(j) * DP(k);
                        A(i,j,k) = (2*i-1) * (2*j-1) * (2*k-1) * int(f, -1, 1);
                    end
                end
            end
            obj.PPDP = A / sqrt(obj.meshSize^3);
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
    methods
        %% deal with DD model
        function [obj, TV, e] = DDModelDGFunction(obj,n)
            x = linspace(0,0.6,10000);
            Time = [];
            TV = [];
            FLast = obj.initialFunction;
            [obj, E, Func] = obj.GK(n+1);
            FNow = Func;
            e = sqrt(gaussLegendre(@(x) (FLast.solve(x) - FNow.solve(x)).^2, obj.xa, obj.xb));
            while  e(end) > obj.epsilon
                tic
                [obj, E, Func] = obj.GK(n+1);
                Time = [Time toc];
                FLast = FNow;
                FNow = Func;
                Y = FNow.solve(obj.Xc);
                TV = [TV sum(abs(Y(2:end)-Y(1:end-1)))];
                e = [e sqrt(gaussLegendre(@(x) (FLast.solve(x) - FNow.solve(x)).^2, obj.xa, obj.xb))];
                %                 obj.save('workspace');
                %                 save('workspace', 'FNow', 'Time', 'TV', '-append');
            end
            disp('over');
        end
        %% GK
        function [obj, E, Func] = GK(obj, n)
            if n == 3
                [obj, E, Func] = obj.GK3;
            elseif n == 4
                [obj, E, Func] = obj.GK4;
            else
                error('Invalid GK order. Please choose 3 or 4.');
            end
        end
        
        function [obj, E, Func] = GK3(obj)
            x = linspace(0,0.6,1000);
            k0 = obj.coeffs;
            F = L(obj);
            k1 = k0 + obj.t * F;
            obj.coeffs = k1;
            F = L(obj);
            k2 = 3/4 * k0 + 1/4 * k1 + 1/4 * obj.t * F;
            obj.coeffs = k2;
            [F,E] = L(obj);
            obj.coeffs = 1/3 * k0 + 2/3 * k2 + 2/3 * obj.t * F;
            Func = LegendrePoly(obj.X,reshape(obj.coeffs,obj.degree+1,[]),obj.degree);
        end
        function [obj,E,Func] = GK4(obj)
            x = linspace(0,0.6,1000);
            k0 = obj.coeffs;
            F0 = L(obj);
            k1 = k0 +  1/ 2 * obj.t * F0;
            obj.coeffs = k1;
            F1 = L(obj);
            k2 = 1/2 * (k0 + k1)  + obj.t*(-1/4*F0 + 1/2*F1);
            obj.coeffs = k2;
            F2 = L(obj);
            k3 = (1*k0 + 2*k1)/9 + 2/3*k2 + obj.t*(-1/9*F0 - 1/3*F1 + F2);
            obj.coeffs = k3;
            [F3,E] = L(obj);
            obj.coeffs = (k1+k2+k3)/3 + obj.t*(F1+F3)/6;
            Func = LegendrePoly(obj.X,reshape(obj.coeffs,obj.degree+1,[]),obj.degree);
        end
        
        % input: obj
        % output: F, electric field E, electron concentration func;
        function [F, E] = L(obj)
            if strcmp(obj.semiModel.name,'DD')
                [F,E] = obj.LDD;
            elseif strcmp(obj.semiModel.name, 'HF')
                [F,E] = obj.LHF;
            else
                error('Invalid model. Please choose either "DD" or "HF".');
            end
        end
        % input mesh info and electron concentration coeffs
        % output electronConcentration, E
        function [electronConcentration, E] = getFuncs(obj)
            electronConcentration = LegendrePoly(obj.X, reshape(obj.coeffs,obj.degree+1,obj.CellsNum),obj.degree,'true');
            temp = electronConcentration.priSolve(0.6);
            T1 = gaussLegendre(@(x) electronConcentration.priSolve(x), 0, 0.6)  + gaussLegendre(@(x) electronConcentration.priSolve(x)+temp,0,0.4);
            temp = obj.initialFunction.priSolve(0.6);
            T2 = gaussLegendre(@(x) obj.initialFunction.priSolve(x), 0, 0.6)  + gaussLegendre(@(x) obj.initialFunction.priSolve(x)+temp,0,0.4); % set in generate function
            electricField0 = 0.001546423010635 * (T1 - T2);
            
            temp = -0.001546423010635 * (-electronConcentration.priSolve(0) + obj.initialFunction.priSolve(0)) + electricField0 -  1.5;
            E = @(x) -0.001546423010635 * (electronConcentration.priSolve(x) - obj.initialFunction.priSolve(x) ) + temp;
        end
        function [F,E] = LDD(obj)
            n = obj.degree;
            N = obj.CellsNum;
            parforX = obj.X;
            parforCoeffs = reshape(obj.coeffs,n+1,N);
            
            temp = sqrt(obj.semiModel.theta);
            sqrtTauTheta = @(x) temp*sqrt(obj.semiModel.relaxationParameter(x));
            

            
%             electronConcentration = LegendrePoly(parforX, parforCoeffs,n,'true');
%             CellValues(:,2:3) = electronConcentration.getNodeValues';
            
            %            Y = electronConcentration.solve(obj.Xc);
            %             tempR = obj.ENO(Y,2,"right");
            %             tempL = obj.ENO(Y,2,"left");
            %             CellValues(:,2) = [tempR(end) tempR(1:end-1)]';
            %             CellValues(:,3) = tempL';
            
            
            % minmod limiter
            M = 2/3 * 1.96e9 * obj.meshSize^2;
            obj.coeffs = obj.minmodLimiter(obj.coeffs, M);
%                         CellAverages = zeros(obj.CellsNum,1);
%                         for j = 1:obj.CellsNum
%                             CellAverages(j) = gaussLegendre(@(x) electronConcentration.solve(x), obj.X(j), obj.X(j+1));
%                         end
%                         CellAverages = CellAverages(:) / obj.meshSize;
%                         M = 2/3 * 2 * 10^7 * obj.meshSize^2;
%                         ur = CellValues(:,3) - CellAverages;
%                         ul = CellAverages - CellValues(:,2);
%                         CellAveragesP = circshift(CellAverages,-1) - CellAverages;
%                         CellAveragesN = CellAverages - circshift(CellAverages,1);
%                         Ar = [ur CellAveragesP CellAveragesN];
%                         Al = [ul CellAveragesP CellAveragesN];
%                         urmod = CellAverages + minmod(Ar,M)';
%                         ulmod = CellAverages - minmod(Al,M)';
%                         index = (abs(urmod-CellValues(:,3)) < 1 ) .* (abs(ulmod - CellValues(:,2)) < 1);
%                         index = ~index;
%                         parforCoeffs(1,index) = sqrt(obj.meshSize) * CellAverages(index)';
%                         parforCoeffs(2,index) = (urmod(index) - ulmod(index)) / (2*sqrt( 3/obj.meshSize ));
%                         parforCoeffs(3,index) = urmod(index) - parforCoeffs(1,index)' / sqrt(obj.meshSize) - parforCoeffs(2,index)' * sqrt(3/obj.meshSize);
            % regenerate
            electronConcentration = LegendrePoly(parforX, reshape(obj.coeffs,n+1,N),n,'true');
            x = linspace(0,0.6,10000);plot(x,electronConcentration.solve(x));
            CellValues = zeros(N,4);%[ll lr rl rr]
            CellValues(:,2:3) = electronConcentration.getNodeValues';
            CellValues(:,1) = circshift(CellValues(:,3),1);
            CellValues(:,4) = circshift(CellValues(:,2),-1);
            
            % get electric field
            temp = electronConcentration.priSolve(0.6);
            T1 = gaussLegendre(@(x) electronConcentration.priSolve(x), 0, 0.6)  + gaussLegendre(@(x) electronConcentration.priSolve(x)+temp,0,0.4);
            temp = obj.initialFunction.priSolve(0.6);
            T2 = gaussLegendre(@(x) obj.initialFunction.priSolve(x), 0, 0.6)  + gaussLegendre(@(x) obj.initialFunction.priSolve(x)+temp,0,0.4); % set in generate function
            electricField0 = 0.001546423010635 * (T1 - T2);
            
            temp = -0.001546423010635 * (-electronConcentration.priSolve(0) + obj.initialFunction.priSolve(0)) + electricField0 -  1.5;
            E = @(x) -0.001546423010635 * (electronConcentration.priSolve(x) - obj.initialFunction.priSolve(x) ) + temp;
            
            % aux function
            
            if isa(obj.semiModel.mobility, 'function_handle')
                sqrtTauThetaValues = [sqrtTauTheta(obj.X(1:end-1)') sqrtTauTheta(obj.X(2:end)')];
                TT1 = getGaussLegendreB(@(x) sqrtTauTheta(x) .* electronConcentration.solve(x), obj.X(1:end-1), obj.X(2:end));
                TT3 = cellfun(@(f) getGaussLegendreB(@(x) f((2*x - obj.X(1:end-1)-obj.X(2:end))/obj.meshSize), obj.X(1:end-1), obj.X(2:end)),obj.diffBasisFuncs, 'UniformOutput', false);
                [~,~,C] = gaussLegendre(@(x) 0*x,obj.X(1:end-1),obj.X(2:end));
                preComp = TT1 .* C;
                result = cellfun(@(x) sum(preComp .* x , 1), TT3, 'UniformOutput', false);
                result = cell2mat(result); %  back to a matrix
                result = reshape(result,[],1);
                
                index = repmat(1:obj.CellsNum,obj.degree+1,1);
                index = reshape(index,[],1);
                F = repmat(obj.basisBoundaryValues,obj.CellsNum,1);
                tempCellValues = CellValues(index,:);
                tempSqrtTauThetaValues = sqrtTauThetaValues(index,:);
                
                obj.auxCoeffs  =tempSqrtTauThetaValues(:,2) .*tempCellValues(:,4) .* F(:,2) -tempSqrtTauThetaValues(:,1).* tempCellValues(:,2) .* F(:,1) - result;
            else
                
                I = eye(obj.CellsNum);
                B = kron(I, obj.PnDPm');
                index = repmat(1:obj.CellsNum,obj.degree+1,1);
                index = reshape(index,[],1);
                
                F = repmat(obj.basisBoundaryValues,obj.CellsNum,1);
                tempCellValues = CellValues(index,:);
                F = tempCellValues(:,4) .* F(:,2) - tempCellValues(:,2) .* F(:,1) - B * obj.coeffs;
                obj.auxCoeffs = 0.139219332249189 * F;
            end
            auxq = LegendrePoly(parforX, reshape(obj.auxCoeffs,obj.degree+1,obj.CellsNum), n);
            x = linspace(0,0.6,10000);plot(x,auxq.solve(x));
            % minimode
            M = 2/3 * 4.4050e+11 * obj.meshSize^2;
            obj.auxCoeffs = obj.minmodLimiter(obj.auxCoeffs, M);
            
            
            auxq = LegendrePoly(parforX, reshape(obj.auxCoeffs,obj.degree+1,obj.CellsNum), n);
            x = linspace(0,0.6,10000);plot(x,auxq.solve(x));
            auxCellValues(:,2:3) = auxq.getNodeValues';
            auxCellValues(:,1) = circshift(auxCellValues(:,3),1);
            auxCellValues(:,4) = circshift(auxCellValues(:,2),-1);
            
            % input: parforCells,CellValues, auxCellValues
            % L output F very slow
            % temp var: Cellj CellValuej auxCellValuej Fj
            
            tbv = repmat(obj.basisBoundaryValues,obj.CellsNum,1);
            
            Eb = E(parforX(2:end));
            Ea = E(parforX(1:end-1));
            if isa(obj.semiModel.mobility, 'function_handle')
                sqrtTauThetaValues = [sqrtTauTheta(obj.X(1:end-1)') sqrtTauTheta(obj.X(2:end)')];
                T3 = obj.semiModel.mobility(obj.X(2:end)) .* (max(Eb,0) .* CellValues(:,4) + min(Eb,0) .* CellValues(:,3)) + sqrtTauThetaValues(:,2) .* auxCellValues(:,3);
                T4 = obj.semiModel.mobility(obj.X(1:end-1)) .* (max(Ea,0) .* CellValues(:,2) + min(Ea,0) .* CellValues(:,1)) + sqrtTauThetaValues(:,1) .* auxCellValues(:,1);
            else
                T3 = 0.75 * (max(Eb,0) .* CellValues(:,4) + min(Eb,0) .* CellValues(:,3)) + 0.139219332249189 * auxCellValues(:,3);
                T4 = 0.75 * (max(Ea,0) .* CellValues(:,2) + min(Ea,0) .* CellValues(:,1)) + 0.139219332249189 * auxCellValues(:,1);
            end
            T3 = T3(index);
            T4 = T4(index);
            
            if isa(obj.semiModel.mobility, 'function_handle')
                TT1 = getGaussLegendreB(@(x) obj.semiModel.mobility(x) .* E(x), obj.X(1:end-1), obj.X(2:end));
            else
                TT1 = 0.75 * getGaussLegendreB(@(x) E(x),obj.X(1:end-1),obj.X(2:end));
            end
            TT2  = getGaussLegendreB(@(x) electronConcentration.solve(x), obj.X(1:end-1), obj.X(2:end));
            TT3 = cellfun(@(f) getGaussLegendreB(@(x) f((2*x - obj.X(1:end-1)-obj.X(2:end))/obj.meshSize), obj.X(1:end-1), obj.X(2:end)),obj.diffBasisFuncs, 'UniformOutput', false);
            [~,~,C] = gaussLegendre(@(x) 0*x,obj.X(1:end-1),obj.X(2:end));
            preComp = TT1 .* TT2 .* C;
            result = cellfun(@(x) sum(preComp .* x , 1), TT3, 'UniformOutput', false);
            result = cell2mat(result); %  back to a matrix
            T1 = reshape(result,[],1);
            if isa(obj.semiModel.mobility, 'function_handle')
                TT1 = getGaussLegendreB(@(x) sqrtTauTheta(x) .* auxq.solve(x), obj.X(1:end-1), obj.X(2:end));
                preComp = TT1 .* C;
                result = cellfun(@(x) sum(preComp .* x , 1), TT3, 'UniformOutput', false);
                result = cell2mat(result); %  back to a matrix
                T2 = reshape(result,[],1);
            else
                T2 = 0.139219332249189 * B * obj.auxCoeffs;
            end
            
            F = T3.*tbv(:,2) - T4.*tbv(:,1) -T2 - T1;
        end
        function f = getLegendrePoly(obj, c)
            f = LegendrePoly(obj.X,reshape(c,obj.degree+1,[]),obj.degree);
        end
        function [F, E] = LHF(obj)
            n = obj.degree;
            N = obj.CellsNum;
            coeffMatrix = reshape(obj.coeffs,n+1,N);
            
            % get electron concentration
            electronConcentration = LegendrePoly(obj.X, coeffMatrix,n,'true');
            
            % get cell values
            CellValues = zeros(N,4);%[ll lr rl rr]
            CellValues(:,2:3) = electronConcentration.getNodeValues';
            CellValues(:,1) = circshift(CellValues(:,3),1);
            CellValues(:,4) = circshift(CellValues(:,2),-1);
            
            % get electric field
            temp = electronConcentration.priSolve(0.6);
            T1 = gaussLegendre(@(x) electronConcentration.priSolve(x), 0, 0.6)  + gaussLegendre(@(x) electronConcentration.priSolve(x)+temp,0,0.4);
            temp = obj.initialFunction.priSolve(0.6);
            T2 = gaussLegendre(@(x) obj.initialFunction.priSolve(x), 0, 0.6)  + gaussLegendre(@(x) obj.initialFunction.priSolve(x)+temp,0,0.4); % set in generate function
            electricField0 = 0.001546423010635 * (T1 - T2);
            temp = -0.001546423010635 * (-electronConcentration.priSolve(0) + obj.initialFunction.priSolve(0)) + electricField0 -  1.5;
            E = @(x) -0.001546423010635 * (electronConcentration.priSolve(x) - obj.initialFunction.priSolve(x) ) + temp;
            
            % get constant
            C = [8.027379234129259e-35 6.020534425596945e-35];
            C(3) = E(0) * electronConcentration.solve(0) * 0.75 * C(1);
            
            % get a b
            a = @(x) sqrt(0.019382022471910 + 1.108773408239700e-32 * 0.75^2 * E(x).^2);
            b = @(x) 3 * C(2) * E(x) .* dopingFunction(x) + 0.75 * E(x) - C(3);
            
            % get auxiliary function
            I = eye(obj.CellsNum);
            B = kron(I, obj.PnDPm');
            index = repmat(1:obj.CellsNum,obj.degree+1,1);
            index = reshape(index,[],1);
            F = repmat(obj.basisBoundaryValues,obj.CellsNum,1);
            tempCellValues = CellValues(index,:);
            [~,~,W] = gaussLegendre(@(x) 0*x,obj.X(1:end-1),obj.X(2:end));
            
            T11 = getGaussLegendreB(@(x) a(x),obj.X(1:end-1),obj.X(2:end));
            T12  = getGaussLegendreB(@(x) electronConcentration.solve(x), obj.X(1:end-1), obj.X(2:end));
            T13 = cellfun(@(f) getGaussLegendreB(@(x) f((2*x - obj.X(1:end-1)-obj.X(2:end))/obj.meshSize), obj.X(1:end-1), obj.X(2:end)),obj.diffBasisFuncs, 'UniformOutput', false);
            preComp = T11 .* T12 .* W;
            result = cellfun(@(x) sum(preComp .* x , 1), T13, 'UniformOutput', false);
            result = cell2mat(result); %  back to a matrix
            T1 = reshape(result,[],1);
            
            dE = @(x) -0.001546423010635 * (electronConcentration.solve(x) - dopingFunction(x));
            da = @(x) (0.019382022471910 + 2*1.108773408239700e-32* 0.75^2 * E(x) .* dE(x)) ./ (2*a(x));
            T21 = getGaussLegendreB(@(x) da(x),obj.X(1:end-1),obj.X(2:end));
            T22 = T12;
            T23 = cellfun(@(f) getGaussLegendreB(@(x) f((2*x - obj.X(1:end-1)-obj.X(2:end))/obj.meshSize), obj.X(1:end-1), obj.X(2:end)),obj.classicalLegendrePolys, 'UniformOutput', false);
            preComp = T21 .* T22 .* W;
            result = cellfun(@(x) sum(preComp .* x , 1), T23, 'UniformOutput', false);
            result = cell2mat(result); %  back to a matrix
            T2 = reshape(result,[],1);
            
            T3 = a(obj.X(2:end)) .* CellValues(:,4) ;
            T4 = a(obj.X(1:end-1)) .* CellValues(:,2) ;
            T3 = T3(index);
            T4 = T4(index);
            obj.auxCoeffs = -T1 - T2 + T3.* F(:,2) - T4.* F(:,1);
            
            auxq = LegendrePoly(obj.X, reshape(obj.auxCoeffs,obj.degree+1,obj.CellsNum), n);
            
            % get auxiliary function Cell Boundary Values
            auxCellValues(:,2:3) = auxq.getNodeValues';
            auxCellValues(:,1) = circshift(auxCellValues(:,3),1);
            auxCellValues(:,4) = circshift(auxCellValues(:,2),-1);
            
            % get electron concentration
            T11 = getGaussLegendreB(@(x) b(x),obj.X(1:end-1),obj.X(2:end));
            preComp = T11 .* T12 .* W;
            result = cellfun(@(x) sum(preComp .* x , 1), T13, 'UniformOutput', false);
            result = cell2mat(result); %  back to a matrix
            T1 = reshape(result,[],1);
            
            idx = b(obj.X(2:end)) >0 ;
            T2 = idx .* b(obj.X(2:end)) .* CellValues(:,4) + ~idx .* b(obj.X(2:end)) .* CellValues(:,3);
            T2 = T2(index);
            T2 = T2  .* F(:,2);
            idx = b(obj.X(1:end-1)) > 0 ;
            T3 = idx .* b(obj.X(1:end-1)) .* CellValues(:,2)  + ~idx .* b(obj.X(1:end-1)) .* CellValues(:,1);
            T3 = T3(index);
            T3 = T3 .* F(:,1);
            
            T41 = gaussLegendre(@(x) 2*C(2) * E(x) .* (electronConcentration.solve(x)).^2,obj.X(1:end-1), obj.X(2:end));
            T43 = T13;
            preComp = T41 .* W;
            result = cellfun(@(x) sum(preComp .* x , 1), T43, 'UniformOutput', false);
            result = cell2mat(result); %  back to a matrix
            T4 = reshape(result,[],1);
            
            T5 = 2 * C(2) * E(obj.X(2:end)) .* (CellValues(:,3).^2 + CellValues(:,3).*CellValues(:,4) + CellValues(:,3).^2 ) ;
            T6 = 2 * C(2) * E(obj.X(2:end)) .* (CellValues(:,1).^2 + CellValues(:,1).*CellValues(:,2) + CellValues(:,1).^2 ) ;
            T5 = T5(index).* F(:,2);
            T6 = T6(index).* F(:,1);
            
            T71 = getGaussLegendreB(@(x) a(x),obj.X(1:end-1),obj.X(2:end));
            T72 = getGaussLegendreB(@(x) auxq.solve(x), obj.X(1:end-1), obj.X(2:end));
            T73 = T13;
            preComp = T71 .* T72 .* W;
            result = cellfun(@(x) sum(preComp .* x , 1), T73, 'UniformOutput', false);
            result = cell2mat(result); %  back to a matrix
            T7 = reshape(result,[],1);
            
            T8 = a(obj.X(2:end)) .* auxCellValues(:,3) ;
            T9 = a(obj.X(1:end-1)) .* auxCellValues(:,1);
            T8 = T8(index).* F(:,2);
            T9 = T9(index) .* F(:,1);
            F = -T1 + T2 - T3 + T4 - T5 + T6 - T7 + T8 - T9;
        end
    end
    methods (Access = private)
        function result = ENO(obj,Y,m,var)
            m = 2*m+1;
            N = obj.CellsNum;
            temp = obj.meshSize * ones(1,m);
            R = cumsum(temp);
            L = cumsum(temp,'reverse');
            X = [obj.Xc(1)-L obj.Xc obj.Xc(end)+R];
            Y = [Y(end-m+1:end) Y Y(1:m)];
            a = zeros(m,N);
            b = zeros(m,N);
            c = zeros(m,N);
            Q = zeros(m,N);
            kmin = zeros(m,N);
            kmax = zeros(m,N);
            if var == "right"
                kmin(1,:) = 1+m:N+m;
                kmax(1,:) = 1+m:N+m;
                Q(1,:) = Y(1+m:N+m);
            elseif var == "left"
                kmin(1,:) = 1+m+1:N+m+1;
                kmax(1,:) = 1+m+1:N+m+1;
                Q(1,:) = Y(1+m+1:N+m+1);
            end
            xj = obj.X(2:1+N); % row vector
            %             table = cell(1,N);
            %             table(:) = {[]};
            for i = 2:m
                % aXs is cloumn cell with aXs(j) being aXj in i state
                aXs = arrayfun(@(a,b) X(a:b), kmin(i-1,:), kmax(i-1,:)+1, 'UniformOutput', false);
                aYs = arrayfun(@(a,b) Y(a:b), kmin(i-1,:), kmax(i-1,:)+1, 'UniformOutput', false);
                aTables = cellfun(@(x,y) dividedDiff(x,y), aXs, aYs, 'UniformOutput',false);
                a(i,:) = cellfun(@(x) x(end,end), aTables);
                
                bXs = arrayfun(@(a,b) X(a:b), kmin(i-1,:)-1, kmax(i-1,:), 'UniformOutput', false);
                bYs = arrayfun(@(a,b) Y(a:b), kmin(i-1,:)-1, kmax(i-1,:), 'UniformOutput', false);
                bTables = cellfun(@(x,y) dividedDiff(x,y), bXs, bYs, 'UniformOutput',false);
                b(i, :) = cellfun(@(x) x(end,end), bTables);
                %                 if i == 2
                %                     % aXs is cloumn cell with aXs(j) being aXj in i state
                %                     aXs = arrayfun(@(a,b) X(a:b), kmin(i-1,:), kmax(i-1,:)+1, 'UniformOutput', false);
                %                     aYs = arrayfun(@(a,b) Y(a:b), kmin(i-1,:), kmax(i-1,:)+1, 'UniformOutput', false);
                %                     aTables = cellfun(@(x,y,z) dividedDiff(x,y,z), aXs, aYs, table, 'UniformOutput',false);
                %                     a(i,:) = cellfun(@(x) x(end,end), aTables);
                %
                %                     bXs = arrayfun(@(a,b) X(a:b), kmin(i-1,:)-1, kmax(i-1,:), 'UniformOutput', false);
                %                     bYs = arrayfun(@(a,b) Y(a:b), kmin(i-1,:)-1, kmax(i-1,:), 'UniformOutput', false);
                %                     bTables = cellfun(@(x,y,table) dividedDiff(x,y,table), bXs, bYs, table, 'UniformOutput',false);
                %                     b(i, :) = cellfun(@(x) x(end,end), bTables);
                %                 else
                %                     aXs = X(anewinterpPoints);
                %                     aYs = Y(anewinterpPoints);
                %                     bXs = X(bnewinterpPoints);
                %                     bYs = Y(bnewinterpPoints);
                %                     aXs = num2cell(aXs);
                %                     aYs = num2cell(aYs);
                %                     bXs = num2cell(bXs);
                %                     bYs = num2cell(bYs);
                %                     aTables = cellfun(@(x,y,z) dividedDiff(x,y,z), aXs, aYs, aTables, 'UniformOutput',false);
                %                     a(i,:) = cellfun(@(x) x(end,end), aTables);
                %                     bTables = cellfun(@(x,y,table) dividedDiff(x,y,table), bXs, bYs, bTables, 'UniformOutput',false);
                %                     b(i, :) = cellfun(@(x) x(end,end), bTables);
                %                 end
                
                temp = abs(a(i,:))+1 >= abs(b(i,:));
                c(i,:) = b(i,:) .* temp + a(i,:) .* (~temp);
                kmin(i,:) = kmin(i-1,:) - temp;
                kmax(i,:) = kmax(i-1,:) + ~temp;
                %                 anewinterpPoints = kmin(i,:) .* temp + (kmax(i,:)+1) .* ~temp;
                %                 bnewinterpPoints = (kmin(i,:)-1) .* temp + kmax(i,:) .* ~temp;
                %                 table = bTables;
                %                 idx = find(~temp);
                %                 for ii=idx
                %                     table{ii} = aTables{ii};
                %                 end
                stencil = arrayfun(@(a,b) X(a:b), kmin(i-1,:), kmax(i-1,:), 'UniformOutput', false);
                stencil = cell2mat(stencil');
                stencil = stencil';
                p = xj - stencil;
                p = prod(p,1);
                
                Q(i,:) = c(i,:) .* p;
            end
            result = sum(Q,1);
        end
        
        
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
        function c = minmodLimiter(obj, c, M)
            c = reshape(c,obj.degree+1,obj.CellsNum);
            f = LegendrePoly(obj.X, c,obj.degree);
            CellValues(:,2:3) = f.getNodeValues';
            CellValues(:,1) = circshift(CellValues(:,3),1);
            CellValues(:,4) = circshift(CellValues(:,2),-1);
            CellAverages = zeros(obj.CellsNum,1);   
            for j = 1:obj.CellsNum
                CellAverages(j) = gaussLegendre(@(x) f.solve(x), obj.X(j), obj.X(j+1));
            end
            CellAverages = CellAverages(:) / obj.meshSize;
            ur = CellValues(:,3) - CellAverages;
            ul = CellAverages - CellValues(:,2);
            CellAveragesP = circshift(CellAverages,-1) - CellAverages;
            CellAveragesN = CellAverages - circshift(CellAverages,1);
            Ar = [ur CellAveragesP CellAveragesN];
            Al = [ul CellAveragesP CellAveragesN];
            urmod = CellAverages + minmod(Ar,M)';
            ulmod = CellAverages - minmod(Al,M)';
%             index = (abs(urmod-CellValues(:,3)) < 1 ) .* (abs(ulmod - CellValues(:,2)) < 1);
%             index = ~index;
            index = 1:obj.CellsNum;
            index = index(:);
            temp = zeros(obj.degree+1,obj.CellsNum);
            temp(1,index) =  CellAverages(index)' * sqrt(obj.meshSize);
            temp(2,index) = (urmod(index) - ulmod(index)) / (2*sqrt( 3/obj.meshSize ));
            temp(3,index) = sqrt(obj.meshSize) * urmod(index) - temp(1,index)' - temp(2,index)' * sqrt(3);
            temp(3,index) = temp(3,index)/sqrt(5);
            c = reshape(temp, [], 1);
        end

    end
    
    
    
end

