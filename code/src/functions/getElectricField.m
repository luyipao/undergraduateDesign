function electricField = getElectricField(priElectronConcentration)
global ELECTRON_CHARGE DIELECTRIC_PERMITTIVITY
coeff = 0.0015; % ELECTRON_CHARGE / DIELECTRIC_PERMITTIVITY;

% E0
T1 = gaussLegendre(@(x) priElectronConcentration(x), 0, 1);
T2 = gaussLegendre(@(x) priDopingFunction(x), 0, 1);
electricField0 = coeff * (T1 - T2);
% E^h
global VOLTAGE_DROP

electricField = @(x) -coeff * (priElectronConcentration(x) - priDopingFunction(x) - priElectronConcentration(0) + priDopingFunction(0)) + electricField0 -  1.5000;
end