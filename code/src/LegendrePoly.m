classdef LegendrePoly
    properties
        X (:,1) double
        coeffs (:,:) double
        C
        basisFuncs
        priBasisFuncs
        normalCoeffs
        priNormalCoeffs
        h (1,1) double
    end
    methods
        function obj = LegendrePoly(mesh, coeffs, degree, C)
            if size(coeffs) ~= [degree+1,length(mesh)-1]
                error('input error: size ');
            end
            obj.X = mesh;
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
            obj.h = abs(mesh(2) - mesh(1));
            normalCoeffs = [sqrt(1/obj.h) sqrt(3/obj.h) sqrt(5/obj.h) sqrt(7/obj.h)];
            obj.priNormalCoeffs = [sqrt(obj.h) sqrt(3*obj.h) sqrt(5*obj.h) sqrt(7*obj.h)] / 2;
            obj.normalCoeffs = normalCoeffs(1:degree+1)';
            obj.basisFuncs = obj.basisFuncs(1:degree+1);
            obj.priBasisFuncs = obj.priBasisFuncs(1:degree+1);
            obj.priNormalCoeffs = obj.priNormalCoeffs(1:degree+1)';
            obj.coeffs = coeffs;
            if nargin == 4
                if ~isrow(C)
                    obj.C = C';
                end
            end
        end
		% output: Cells value in both side, y(1,j)=Celljfunc(a), y(2,j) = Celljfunc(b) 
		function y = getNodeValues(obj)
			normalizeCoeffs = [obj.normalCoeffs; obj.normalCoeffs .*(-1).^(0:degree)]
			normalizeCoeffs = repmat(normalizeCoeffs,1,length(obj.X));
			coeffs = repmat(obj.coeffs,2,1);
			A = coeffs .* normalizeCoeffs;
			y(1,:) = sum(A(1:obj.degree+1,:),1);
			y(2,:) = sum(A(obj.degree+2:end,:),1);
		end
        function y = solve(obj, x)
            if isrow(x)
                x = x';
            end
            index = discretize(x,obj.X);
            normalizeCoeffs = repmat(obj.normalCoeffs,1,length(x));
            y = obj.coeffs(:,index) .* normalizeCoeffs .* cell2mat(cellfun(@(f) f((2*x-obj.X(index)-obj.X(index+1)) / (obj.h)), obj.priBasisFuncs, 'UniformOutput', false))';
            y = sum(y,1)';
        end
        function y = priSolve(obj,x)
            if isrow(x)
                x = x';
            end
            index = discretize(x,obj.X);
            priNormalizeCoeffs = repmat(obj.priNormalCoeffs,1,length(x));
            y = obj.coeffs(:,index) .* priNormalizeCoeffs .* cell2mat(cellfun(@(f) f((2*x-obj.X(index)-obj.X(index+1)) / (obj.h)), obj.priBasisFuncs, 'UniformOutput', false))';
            y(1,:) = y(1,:) + obj.C(index);
            y = sum(y,1)';
            % add constant value
        end
    end
end