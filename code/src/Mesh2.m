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
        BPos
        BNeg
        
        assMatrix
        
        CFL
        epsilon = 0.001
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
            
            obj = obj.setTimeStep(1.2e-3);
            obj = obj.setBNegBpos;
            obj = obj.setAssMatrix;
            
            obj = obj.L2projection(f);
            obj = obj.setInitialCoeffs;
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
                f3 = @(x) x.^3 - 3/20 * x;
                obj.basisFuncs = {f0 f1 f2 f3};
                obj.basisFuncs = obj.basisFuncs(1:obj.degree+1);
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
        function obj = setBNegBpos(obj)
            obj.BPos = obj.getHPos;
            obj.BNeg = obj.getHNeg;
        end
        function obj = setAssMatrix(obj)
            A = obj.meshSize * kron(eye(obj.CellsNum),obj.massMatrix);
            obj.assMatrix = A - obj.t/2 * obj.BNeg*diag(1./diag(A))*obj.BPos;
        end
        function obj = setTimeStep(obj, timeStep)
            obj.t = timeStep;
        end
        function f = getBasisPolys(obj,c)
            f = basisPolys(obj.X, reshape(c,obj.degree+1,obj.CellsNum), obj.degree, obj.basisFuncs, obj.basisInterval);
        end
        function [obj, nSteps] = IMEXGK(obj, n)
            
            if n ~= 3
                error('Invalid input, please choose 3');
            else
                % 1
                M = obj.meshSize * kron(eye(obj.CellsNum),obj.massMatrix);
                temp = diag(1./diag(M));
                BNegBPos = obj.BNeg*temp*obj.BPos;
                FLast = obj.getBasisPolys(zeros(obj.CellsNum*(obj.degree+1),1));
                FNow = obj.getBasisPolys(obj.coeffs);
                nSteps = 0;
                filename1 = sprintf('DDIMEXRK3Degree3meshCells%dElectronConcentration_t%g.gif', obj.CellsNum, obj.t);
                filename2 = sprintf('DDIMEXRK3Degree3meshCells%d_E_t%g.gif', obj.CellsNum, obj.t);
                fig1 = figure('Visible', 'off');
                fig2 = figure('Visible', 'off');
                while gaussLegendre(@(x) (FLast.solve(x) - FNow.solve(x)).^2, obj.X(1), obj.X(end)) > (obj.epsilon)^2
                    obj = L(obj, M, BNegBPos);
                    FLast = FNow;
                    FNow = obj.getBasisPolys(obj.coeffs);
                    nSteps = nSteps + 1;
                    
                    % save the change of electron concentration and electric field gif
                    figure(fig1)
                    x = linspace(0,0.6,1000);
                    electronConcentration = obj.getBasisPolys(obj.coeffs);
                    plot(x,electronConcentration.solve(x));
                    imind1 = frame2im(getframe(gcf));
                    [imind1, colormap1] = rgb2ind(imind1, 256);
                    if nSteps == 1
                        imwrite(imind1, colormap1,  fullfile('..\docs\images',filename1), 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
                    else
                        imwrite(imind1, colormap1, fullfile('..\docs\images',filename1), 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
                    end
                    
                    figure(fig2)
                    E = obj.getBasisPolys(obj.Ecoeffs);
                    plot(x,E.solve(x));
                    imind2 = frame2im(getframe(gcf));
                    [imind2, colormap2] = rgb2ind(imind2, 256);
                    if nSteps == 1
                        imwrite(imind2, colormap2,  fullfile('..\docs\images',filename2), 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
                    else
                        imwrite(imind2, colormap2, fullfile('..\docs\images',filename2), 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
                    end
                end
            end
        end
        function obj = L(obj, M, BNegBPos)
            H = cell(4,1);
            c = cell(4,1);
            c{1} = obj.coeffs;
            %                 electronConcentration = obj.getBasisPolys(obj.coeffs);
            %                 x = linspace(0,0.6,1000);plot(x,dopingFunction(x), x, electronConcentration.solve(x),'--');
            
            obj = obj.getIMEXPotential;
            H{1} = obj.H;
            obj.coeffs = obj.assMatrix \ (M*c{1} + obj.t/2*H{1});
            c{2} = obj.coeffs;
            
            obj = obj.getIMEXPotential;
            H{2} =  obj.H;
            b = (M*c{1} +(11*H{1}+H{2})*obj.t/18 + BNegBPos*(obj.t/6*c{2}));
            obj.coeffs = obj.assMatrix \ b;
            c{3} = obj.coeffs;
            
            obj = obj.getIMEXPotential;
            H{3} = obj.H;
            b = M*c{1}+obj.t*((5*H{1}-5*H{2}+3*H{3})/6 + BNegBPos*(-c{2}+c{3})/2);
            obj.coeffs = obj.assMatrix \ b;
            c{4} = obj.coeffs;
            
            obj = obj.getIMEXPotential;
            H{4} = obj.H;
            b = (M*c{1}+obj.t*((H{1}+7*H{2}+3*H{3}-7*H{4})/4 + BNegBPos*(3*(c{2}-c{3})+c{4})/2));
            obj.coeffs = obj.assMatrix \ b;
            obj.auxCoeffs = obj.BPos*obj.coeffs;
        end
        function [obj, CellValues, ECellValues] = getCellValues(obj)
            electronConcentration = obj.getBasisPolys(obj.coeffs);
            % get cell values
            CellValues(:,2:3) = electronConcentration.getNodeValues';
            CellValues(:,1) = circshift(CellValues(:,3),1);
            CellValues(:,4) = circshift(CellValues(:,2),-1);
            
            [obj.Ecoeffs, ~] = obj.getIMEXPotential;
            E = obj.getBasisPolys(obj.Ecoeffs);
            ECellValues(:,2:3) = E.getNodeValues';
            ECellValues(:,1) = circshift(ECellValues(:,3),1);
            ECellValues(:,4) = circshift(ECellValues(:,2),-1);
        end
        % output:
        function y = H(obj)
            electronConcentration = obj.getBasisPolys(obj.coeffs);
            nValues(:,2:3) = electronConcentration.getNodeValues';
            nValues(:,1) = circshift(nValues(:,3),1);
            nValues(:,4) = circshift(nValues(:,2),-1);
            
            E = obj.getBasisPolys(obj.Ecoeffs);
            EValues(:,2:3) = E.getNodeValues';
            EValues(:,1) = circshift(EValues(:,3),1);
            EValues(:,4) = circshift(EValues(:,2),-1);
            
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
                    A = obj.PPDP(:,:,l) .* (obj.Ecoeffs(x) * obj.coeffs(x)');
                    result((j-1)*(1+obj.degree)+l) = sum(A,'all');
                end
            end
            y = -result + FEnv1 - FEnv2;
            y = 0.75 * y;
        end
        % B only depend on basis functions values at each Cells.
        function A = getHPos(obj)
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
        
        function A = getHNeg(obj)
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
        
        function [obj, z, x] = getIMEXPotential(obj)
            % get relationship
            [A, b] = obj.IMEXPossionRelation;
            
            % get electric field z and electric potential  x
            I = eye(obj.CellsNum);
            I = sparse(I);
            M = obj.meshSize * kron(I, obj.massMatrix);
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
            
            temp1 = EA * diag(1./diag(M)) * A + PA;
%             electronConcentration = obj.getBasisPolys(obj.coeffs);
%             temp = cellfun(@(r) gaussLegendre(@(x) (electronConcentration.solve(x) - dopingFunction(x)) .* r((x - obj.Xc(1:end))/obj.meshSize),obj.X(1:end-1), obj.X(2:end)), obj.basisFuncs,'UniformOutput', false);
%             temp = cell2mat(temp');
%             temp = reshape(temp,[],1);
%             temp2 = -Pb - EA*(M \ b) - 0.001546423010635 * temp ;
            temp2 = -Pb - EA*(M \ b) - 0.001546423010635 * M * (obj.coeffs- obj.initialCoeffs);
            x =  temp1 \ temp2;
            z =  M \ (A*x + b);
            obj.Ecoeffs = z;
        end
    end
end