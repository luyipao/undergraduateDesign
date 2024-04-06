classdef Cell
    properties
        a (1,1) double
        b (1,1) double
        ll (1,1) double
        lr (1,1) double
        rl (1,1) double
        rr (1,1) double
        auxll (1,1) double
        auxlr (1,1) double
        auxrl (1,1) double
        auxrr (1,1) double
        auxFunction
        coeffs (:,1) double
        auxCoeffs (:,1) double
        % first colum Pi  second  DPi
        basisFunctions (:,2) cell
        basisFunctionsBoundaryValues (:,2) double
        degree double
    end
    methods (Access = public)
        function obj = Cell(a,b,n,coeffs)
            if nargin == 0
                ...
            elseif nargin == 3
            obj.a = a;
            obj.b = b;
            obj.degree = n;
            for i = 1:n+1
                [Pi,DPi] = obj.getBasisFunction(i-1);
                obj.basisFunctions(i,:) = {Pi, DPi};
               obj.basisFunctionsBoundaryValues(i,:) = [Pi(a) Pi(b)];
            end
            elseif nargin == 4
                obj.a = a;
                obj.b = b;
                obj.degree = n;
                obj.coeffs = coeffs;
                for i = 1:n+1
                    [Pi,DPi] = obj.getBasisFunction(i-1);
                    obj.basisFunctions(i,:) = {Pi, DPi};
                    obj.basisFunctionsBoundaryValues(i,:) = [Pi(a) Pi(b)];
                end
            else
                error('Invalid number of inputs.')
            end
        end
    end
    methods (Access = public)
        function [PN,DPN] = getBasisFunction(obj,n)
            P0 = @(x) 0.*x + 1;
            DP0 = @(x) 0.*x + 0;
            P1 = @(x) 1.*x;
            DP1 = @(x) 0.*x + 1;
            if n==0
                P = P0;
                DP = DP0;
            elseif n == 1
                P = P1;
                DP = DP1;
            else
                for ii = 2:n
                    P = @(x) ((2*ii - 1) .* x .* P1(x) - (ii - 1) .* P0(x)) / ii;
                    DP = @(x) ((2*ii - 1) .* P1(x) + (2*ii - 1) .* x .* DP1(x) - (ii - 1) .* DP0(x)) / ii;
                    P0 = P1;
                    DP0 = DP1;
                    P1 = P;
                    DP1 = DP;
                end
            end
            PN = @(x) sqrt((2*n + 1) / (obj.b - obj.a)) * P((2.*x - obj.b - obj.a) / (obj.b - obj.a));
            DPN = @(x) sqrt(4*(2*n + 1) / (obj.b - obj.a)^3) * DP((2.*x - obj.b - obj.a) / (obj.b - obj.a));
        end
    end
    
end