classdef LegendrePoly
    properties
        X (:,1) double
        coeffs (:,:) double
        degree (1,1) 
        C
        diffBasisFuncs
        basisFuncs
        priBasisFuncs
        diffNormalCoeffs
        normalCoeffs
        priNormalCoeffs
        h (1,1) double
        nodeValues
    end
    methods
        function obj = LegendrePoly(mesh, coeffs, degree, ~)
            if size(coeffs) ~= [degree+1,length(mesh)-1]
                error('input error: size ');
            end
            obj.X = mesh;
            obj.degree = degree;
            f0 = @(x) 0*x + 1;
            f1 = @(x) x;
            f2 = @(x) 1.5 * x.^2 - 0.5;
            f3 = @(x) 2.5 * x.^3 - 1.5 * x;
            obj.basisFuncs = {f0 f1 f2 f3};
            f0 = @(x) x;
            f1 = @(x) 1/2 * x.^2;
            f2 = @(x) 1/2 * (x.^3 - x);
            f3 = @(x) 5/8 * x.^4 - 3/4 * x.^2;
            obj.priBasisFuncs = {f0 f1 f2 f3};
            f0 = @(x) 0*x;
            f1 = @(x) 0*x + 1;
            f2 = @(x) 3*x;
            f3 = @(x) 7.5*x - 1.5;
            obj.diffBasisFuncs = {f0 f1 f2 f3};
            
            obj.h = abs(mesh(2) - mesh(1));
            normalCoeffs = [sqrt(1/obj.h) sqrt(3/obj.h) sqrt(5/obj.h) sqrt(7/obj.h)];
            
            obj.diffNormalCoeffs = 2 * normalCoeffs .* [1/obj.h 3/obj.h 5/obj.h 7/obj.h ];
            obj.priNormalCoeffs = [sqrt(obj.h) sqrt(3*obj.h) sqrt(5*obj.h) sqrt(7*obj.h)] / 2;
                        obj.diffNormalCoeffs = obj.diffNormalCoeffs(1:degree+1)';
                        obj.normalCoeffs = normalCoeffs(1:degree+1)';
                        obj.priNormalCoeffs = obj.priNormalCoeffs(1:degree+1)';
            obj.basisFuncs = obj.basisFuncs(1:degree+1);
            obj.priBasisFuncs = obj.priBasisFuncs(1:degree+1);

            obj.coeffs = coeffs;
            obj.nodeValues = obj.getNodeValues;
            if nargin == 4
                temp = obj.getPriNodeValues;
                temp1 = temp(2,:) - temp(1,:);
                temp1 =  [0 temp1];
                temp1 = cumsum(temp1);
                obj.C = temp1(1:end-1) - temp(1,:);
            end
        end
		% output: Cells value in both side, y(1,j)=Celljfunc(a), y(2,j) = Celljfunc(b) 
		function y = getNodeValues(obj)
			normalizeCoeffs = [obj.normalCoeffs; obj.normalCoeffs .*(-1).^(0:obj.degree)'];
			normalizeCoeffs = repmat(normalizeCoeffs,1,length(obj.X)-1);
			A = repmat(obj.coeffs,2,1) .* normalizeCoeffs;
            y(1,:) = sum(A(obj.degree+2:end,:),1);
			y(2,:) = sum(A(1:obj.degree+1,:),1);
        end
		% output cell average value
		function getCellAverages(obj)
			
		end
        % output: y(1,:) left values; y(2,:) right vlaues
		function y = getPriNodeValues(obj)
            temp = [-1 1
                0.5 0.5
                0 0
                -17/24 -17/24];
            temp = temp(1:obj.degree+1,:);
			normalizeCoeffs = [obj.priNormalCoeffs .* temp(:,1); obj.priNormalCoeffs .* temp(:,2)];
			normalizeCoeffs = repmat(normalizeCoeffs,1,length(obj.X)-1);
			A = repmat(obj.coeffs,2,1) .* normalizeCoeffs;
			y(1,:) = sum(A(1:obj.degree+1,:),1);
			y(2,:) = sum(A(obj.degree+2:end,:),1);
        end
        function y = diffSolve(obj, x)
            [m,n] = size(x);
            x = x(:);
            index = discretize(x,obj.X);
            diffNormalizeCoeffs = repmat(obj.diffNormalCoeffs,1,length(x));
            y = obj.coeffs(:,index) .* diffNormalizeCoeffs .* cell2mat(cellfun(@(f) f((2*x-obj.X(index)-obj.X(index+1)) / (obj.h)), obj.diffBasisFuncs, 'UniformOutput', false))'; 
            y = reshape(sum(y,1), m, n);
        end
        function y = solve(obj, x)
            [m,n] = size(x);
            x = x(:);
            
            index = discretize(x,obj.X);
            normalizeCoeffs = repmat(obj.normalCoeffs,1,length(x));
            y = obj.coeffs(:,index) .* normalizeCoeffs .* cell2mat(cellfun(@(f) f((2*x-obj.X(index)-obj.X(index+1)) / (obj.h)), obj.basisFuncs, 'UniformOutput', false))';
            
            y = reshape(sum(y,1), m, n);
        end
        % output y  is row vector
        function y = priSolve(obj,x)
            [m,n] = size(x);
            x = x(:);
            
            index = discretize(x,obj.X);
            priNormalizeCoeffs = repmat(obj.priNormalCoeffs,1,length(x));
            y = obj.coeffs(:,index) .* priNormalizeCoeffs .* cell2mat(cellfun(@(f) f((2*x-obj.X(index)-obj.X(index+1)) / (obj.h)), obj.priBasisFuncs, 'UniformOutput', false))';
            y(1,:) = y(1,:) + obj.C(index);
            
            y = reshape(sum(y,1), m, n);
        end
    end
end