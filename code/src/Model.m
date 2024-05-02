classdef Model
    properties
        name string
        electronCharge (1,1) double
        latticeTemperature (1,1) double
        dielectricPermittivity (1,1) double
        electronEffectiveMass (1,1) double
        BoltzmannConstant (1,1) double
        dopingFunction function_handle
        theta (1,1) double
        mobility
        relaxationParameter
    end
    methods
        function obj = Model(jsonfile)
            data = jsondecode(fileread(jsonfile));
			obj.name = data.name;
            obj.electronCharge = data.electronCharge;
            obj.electronEffectiveMass = data.electronEffectiveMass;
			obj.dielectricPermittivity = data.dielectricPermittivity;
			obj.BoltzmannConstant = data.BoltzmannConstant;
            obj.latticeTemperature = data.latticeTemperature;
            obj.dopingFunction = str2func(data.dopingFunction);
            
            if isa(data.mobility, 'double')
                obj.mobility = data.mobility;
                obj.relaxationParameter = obj.electronEffectiveMass/obj.electronCharge * obj.mobility;
            else
                obj.mobility = str2func(data.mobility);
                temp = obj.electronEffectiveMass/obj.electronCharge;
                obj.relaxationParameter = @(x) temp * obj.mobility(x);
            end
            obj.theta = obj.BoltzmannConstant*obj.latticeTemperature / obj.electronEffectiveMass;
        end
    end
end