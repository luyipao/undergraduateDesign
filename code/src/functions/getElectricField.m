function electricField = getElectricField(priElectronConcentration)
global ELECTRON_CHARGE DIELECTRIC_PERMITTIVITY
coeff = ELECTRON_CHARGE / DIELECTRIC_PERMITTIVITY;

% E0
T1 = quadgk(@(x) priElectronConcentration(x), 0, 1);
T2 = quadgk(@(x) priDopingFunction(x), 0, 1);
electricField0 = coeff * (T1 - T2);
% E^h
global VOLTAGE_DROP

electricField = @(x) -coeff * (priElectronConcentration(x) - priDopingFunction(x) - priElectronConcentration(0) + priDopingFunction(0)) + electricField0 - VOLTAGE_DROP;
end