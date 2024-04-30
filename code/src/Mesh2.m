classdef Mesh2
    properties
        degree (1,1)
        meshSize (1,1)
        CellsNum (1,1)
        X (1,:)
        Xc (1,:)
        coeffs (:,1)
        auxCoeffs (:,1) double
        Ecoeffs (:,1) double
        initialCoeffs (:,1) double
        
        basisFuncs
        basisInterval (1,2) double
        basisBoundaryValues
        massMatrix
        PDP
        PPDP
        
        CFL
        epsilon
        t
    end
    methods
        function obj = Mesh2(a, b, N, f, degree)
            obj.meshSize = (b-a) / N;
            obj.CellsNum = N;
            obj.degree = degree;
            obj.X = linspace(a, b, N+1);
            obj.Xc = (obj.X(1:end-1) + obj.X(2:end)) / 2;
            
            obj = obj.setBasisFuncs;
            obj = obj.setBasisBoundaryValues;
            obj = obj.setMassMatrix;
            obj = obj.setPDP;
            obj = obj.setPPDP;
            
            obj = obj.L2projection(f);
            obj = obj.setInitialCoeffs;
            f = obj.getBasisPolys(obj.coeffs);
            x = linspace(0,0.6,1000);plot(x,f.solve(x));
            

        end

        function obj = setInitialCoeffs(obj)
            c = obj.coeffs;
            obj.initialCoeffs = c;
        end
        function obj = setBasisFuncs(obj, funcs)
            if nargin == 2
                obj.basisFuncs = funcs;
            else
                f0 = @(x) 1+0*x;
                f1 = @(x) x;
                f2 = @(x) x.^2 - 1/12;
                obj.basisFuncs = {f0 f1 f2};
                obj.basisInterval = [-0.5, 0.5];
            end
        end
        % input: projected function f; mesh; basis functions
        % output: f's coeffs on basis functions
        function obj = L2projection(obj, f)
            
            I = eye(obj.CellsNum);
            I = sparse(I);
            A = kron(I, obj.massMatrix);
            
            b = zeros(obj.CellsNum*(obj.degree+1),1);
            for j = 1:obj.CellsNum
                for i = 1:obj.degree+1
                    idx = (j-1)*(obj.degree+1)+i;
                    b(idx) = gaussLegendre(@(x) f(x) .* obj.basisFuncs{i}((x-obj.Xc(j))/obj.meshSize), obj.X(j), obj.X(j+1));
                end
            end
            obj.coeffs = (obj.meshSize * A) \ b;
        end
        function obj = setPDP(obj)
            syms x
            P = sym(obj.basisFuncs);
            DP = diff(P);
            obj.PDP = zeros(obj.degree+1, obj.degree+1);
            for i = 1:obj.degree+1
                for j = 1:obj.degree+1
                    f = P(i) * DP(j);
                    obj.PDP(i,j) = int(f,obj.basisInterval(1),obj.basisInterval(2));
                end
            end
            obj.PDP = obj.PDP;
        end
        function obj = setPPDP(obj)
            syms x
            A = zeros(obj.degree+1, obj.degree+1, obj.degree+1);
            P = sym(obj.basisFuncs);
            DP = diff(P);
            for i = 1:obj.degree+1
                for j = 1:obj.degree+1
                    for k = 1:obj.degree+1
                        f = P(i) * P(j) * DP(k);
                        A(i,j,k) = int(f, obj.basisInterval(1),  obj.basisInterval(2));
                    end
                end
            end
            obj.PPDP = A;
        end
        function obj = setMassMatrix(obj)
            A = zeros(obj.degree+1);
            for i = 1:obj.degree+1
                for j = i
                    g = @(x) obj.basisFuncs{i}(x).*obj.basisFuncs{j}(x);
                    A(i,j) = gaussLegendre(@(x) g(x), obj.basisInterval(1), obj.basisInterval(2));
                end
            end
            obj.massMatrix = A;
        end
        function obj = setBasisBoundaryValues(obj)
            obj.basisBoundaryValues = zeros(obj.degree+1, 2);
            for j = 1:obj.degree+1
                obj.basisBoundaryValues(j,:) = obj.basisFuncs{j}(obj.basisInterval);
            end
        end
        
        function f = getBasisPolys(obj,c)
            f = basisPolys(obj.X, reshape(c,obj.degree+1,obj.CellsNum), obj.degree, obj.basisFuncs, obj.basisInterval);
        end
        function IMEXGK(obj, n)
            obj.t = 1.2e-3;
            if n ~= 3
                error('Invalid input, please choose 3');
            else
                electronConcentration = obj.getBasisPolys(obj.coeffs);
                % get cell values
                CellValues(:,2:3) = electronConcentration.getNodeValues';
                CellValues(:,1) = circshift(CellValues(:,3),1);
                CellValues(:,4) = circshift(CellValues(:,2),-1);
                
                [obj.Ecoeffs, ~] = obj.getIMEXPotential;
                E = basisPolys(obj.X, reshape(obj.Ecoeffs,obj.degree+1,obj.CellsNum), obj.degree, obj.basisFuncs, obj.basisInterval);
                x = linspace(0,0.6,10000); plot(x,E.solve(x))
                ECellValues(:,2:3) = E.getNodeValues';
                ECellValues(:,1) = circshift(ECellValues(:,3),1);
                ECellValues(:,4) = circshift(ECellValues(:,2),-1);
                
                I = sparse(eye(obj.CellsNum*(obj.degree+1)));
                BPos = obj.Hpos;
                BNeg = obj.Hneg;
                
                massMatrix = kron(eye(obj.CellsNum),obj.massMatrix);
                
                H = cell(4,1);
                HPos = cell(4,1);
                HNeg = cell(4,1);
                c = cell(4,1);
                c{1} = obj.coeffs;
                H{1} = obj.H(obj.Ecoeffs, obj.coeffs, ECellValues, CellValues);

                c1 = ((obj.meshSize * massMatrix) -obj.t * BNeg*BPos) \ ((obj.meshSize * massMatrix) * c{1} + obj.t * H{1});
                
                electronConcentration = obj.getBasisPolys(c1);
                x = linspace(0,0.6,1000);plot(x,electronConcentration.solve(x));
            end
            
        end
        % output:
        function y = H(obj, Ecoeffs, ncoeffs, EValues, nValues)
            rbbv = repmat(obj.basisBoundaryValues, obj.CellsNum, 1);
            
            index = repmat(1:obj.CellsNum, obj.degree+1, 1);
            index = reshape(index, [], 1);
            
            EnFlux1 = 0.5 * (EValues(:,4) .* nValues(:,4) + EValues(:,3) .* nValues(:,3));
            FEnv1 = EnFlux1(index) .* rbbv(:,2);
            
            EnFlux2 = 0.5 * (EValues(:,2).*nValues(:,2) + EValues(:,1).*nValues(:,1));
            FEnv2 = EnFlux2(index) .* rbbv(:,1);
            
            result = zeros(obj.CellsNum*(1+obj.degree),1);
            for j = 1:obj.CellsNum
                for l = 1:obj.degree+1
                    x = (j-1)*(1+obj.degree)+1:(j-1)*(1+obj.degree)+obj.degree+1;
                    A = obj.PPDP(:,:,l) .* (Ecoeffs(x) * ncoeffs(x)');
                    result((j-1)*(1+obj.degree)+l) = sum(A,'all');
                end
            end
            
            y = -result + FEnv1 - FEnv2;
            y = 0.75 * y;
        end
        % B only depend on basis functions values at each Cells.
        function A = Hpos(obj)
            I = eye(obj.CellsNum);
            I = sparse(I);
            B = kron(I, obj.PDP');
            
            D = obj.basisBoundaryValues(:,2) * obj.basisBoundaryValues(:,1)';
            D = kron(I, D);
            D = circshift(D,-obj.degree-1,1);
            
            E = obj.basisBoundaryValues(:,1) * obj.basisBoundaryValues(:,1)';
            E = kron(I, E);
            A = -B + D - E;
            A = 0.139219332249189 * A;
        end
        
        function A = Hneg(obj)
            I = eye(obj.CellsNum);
            I = sparse(I);
            B = kron(I, obj.PDP');
            
            D = obj.basisBoundaryValues(:,2) * obj.basisBoundaryValues(:,2)';
            D = kron(I, D);
            
            E = obj.basisBoundaryValues(:,1) * obj.basisBoundaryValues(:,2)';
            E = kron(I, E);
            E = circshift(E, obj.degree+1,1);
            
            A = -B + D - E;
            A = 0.139219332249189 * A;
        end
        function [A, b] = IMEXPossionRelation(obj)
            % get relationship bettween field and potential obj.meshSize * massMatrix * obj.Ecoeffs= A * obj.Pcoeffs + b
            I = eye(obj.CellsNum);
            I = sparse(I);
            B = kron(I, obj.PDP');
            
            temp = obj.basisBoundaryValues(:,2) * obj.basisBoundaryValues(:,1)';
            I(1,1) = 0;
            D = kron(I,temp);
            D = circshift(D, -obj.degree-1, 1);
            
            E = obj.basisBoundaryValues(:,1) * obj.basisBoundaryValues(:,1)';
            I(1,1) = 0;
            E = kron(I,E);
            
            temp = zeros(obj.CellsNum,1);
            temp = sparse(temp);
            temp(end) = 1;
            b = -kron(temp,1.5*obj.basisBoundaryValues(:,2));
            
            A = B - D + E;
        end
        
        function [z, x] = getIMEXPotential(obj)
            % get relationship
            [A, b] = obj.IMEXPossionRelation;
            
            % get electric field z and electric potential  x
            I = eye(obj.CellsNum);
            I = sparse(I);
            massMatrix = kron(I, obj.massMatrix);
            D = obj.basisBoundaryValues(:,2) * obj.basisBoundaryValues(:,2)';
            G = obj.basisBoundaryValues(:,2) * obj.basisBoundaryValues(:,1)';
            F = obj.basisBoundaryValues(:,1) * obj.basisBoundaryValues(:,1)';
            
            AE = kron(I,-obj.PDP'+D);
            BE = kron(diag(sparse(ones(obj.CellsNum-1,1)), -1), G');
            BE(1:obj.degree+1, 1:obj.degree+1) = F;
            EA = AE - BE;
            
            AP = kron(sparse(diag(ones(obj.CellsNum-1,1), 1)),G);
            BP = kron(I,-D-F);
            CP = kron(sparse(diag(ones(obj.CellsNum-1,1), -1)),G');
            PA = AP + BP + CP;
            
            Pb = zeros(obj.CellsNum*(obj.degree+1),1);
            Pb(end-obj.degree:end) = 1.5 * obj.basisBoundaryValues(:,2);
            
            temp1 = EA * inv(obj.meshSize * massMatrix) * A + PA;
            
            %electronConcentration = LegendrePoly(obj.X, reshape(obj.coeffs,obj.degree+1,obj.CellsNum), obj.degree);
            temp = cellfun(@(r) gaussLegendre(@(x) dopingFunction(x) .* r((x - obj.Xc(1:end))/obj.meshSize),obj.X(1:end-1), obj.X(2:end)), obj.basisFuncs,'UniformOutput', false);
            temp = cell2mat(temp');
            temp = reshape(temp,[],1);
            temp = (obj.meshSize * massMatrix) \ temp;
            temp2 = -Pb - EA*(obj.meshSize * massMatrix \ b) - 0.001546423010635 * obj.meshSize * massMatrix * (obj.coeffs- temp);
            
            x =  (temp1 \ temp2);
            z =  (obj.meshSize * massMatrix) \ (A*x + b);
            potential = basisPolys(obj.X, reshape(x,obj.degree+1,obj.CellsNum), obj.degree, obj.basisFuncs, obj.basisInterval);
            x = linspace(0,0.6,10000); plot(x,potential.solve(x))
        end
    end
end