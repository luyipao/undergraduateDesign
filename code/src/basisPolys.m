classdef basisPolys
    properties
        X (:,1) double
        Xc (:,1) double
        h
        coeffs (:,:) double
        degree (1,1)
        basisFuncs
        basisInterval
		basisBoundaryValues
    end
    methods
        function obj = basisPolys(mesh, coeffs, degree, basisFuncs, interval)
            if size(coeffs) ~= [degree+1,length(mesh)-1]
                error('input error: size ');
            else
                obj.X = mesh;
                obj.Xc = 0.5*(mesh(1:end-1) + mesh(2:end));
                obj.degree = degree;
                obj.basisFuncs = basisFuncs;
                obj.basisInterval = interval;
				obj.basisBoundaryValues = zeros(obj.degree+1, 2);
				for j = 1:obj.degree+1
					obj.basisBoundaryValues(j,:) = obj.basisFuncs{j}(obj.basisInterval);
				end
				obj.h = abs(mesh(2) - mesh(1));
                
                obj.coeffs = coeffs;
            end
        end
        function y = getNodeValues(obj)
			A = repmat(obj.coeffs,2,1).*obj.basisBoundaryValues(:);
			y(1,:) = sum(A(1:obj.degree+1,:),1);
            y(2,:) = sum(A(obj.degree+2:end,:),1);
			
		end
        function y = solve(obj, x)
            [m,n] = size(x);
            x = x(:);
			% scale
			[x, index] = transform(obj, x);
			
            y = obj.coeffs(:,index) .* cell2mat(cellfun(@(f) f(x), obj.basisFuncs, 'UniformOutput', false))';
            
            y = reshape(sum(y,1), m, n);
        end

        % transform function: Projects x onto the equivalent points on the support interval of the basis functions
        % Input: x is the vector of points to be transformed
        % Output: x is the vector of transformed points
        function [x, index] = transform(obj, x)
            index = discretize(x,obj.X);
            x = (x-obj.Xc(index))/obj.h;
        end
    end
    
end