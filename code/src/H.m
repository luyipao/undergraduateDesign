% input mobility; electirc field; electron concentration, test function, CellValues, CellElectricFieldValues
function y = H(obj, E, n, v, CellElectricFieldValues, CellValues)
T1 = gaussLegendre(mu .* E(x) .* n(x) .* v(x), obj.X(1:end-1), obj.X(2:end));
T2 = 0.5 * mu .* (CellElectricFieldValues(:,3).*CellValues(:,3) + CellElectricFieldValues(:,4).*CellValues(:,4)) .* BBV(:,2);
T3 = 0.5 * mu .* (CellElectricFieldValues(:,1).*CellValues(:,1) + CellElectricFieldValues(:,2).*CellValues(:,2)) .* BBV(:,1);
y = -T1 + T2 - T3;
end

function y = Hl(obj, u, v, uValues, vValues, dir)
T1 = -sqrt(tau*theta) * gaussLegendre(@(x) u(x) .* v(x), obj.X(1:end-1), obj.X(2:end));
if strcmp(dir, 'pos')
	T2 = sqrt(tau*theta) * uValues(:,4) .* vValues(:,2);
	T3 = sqrt(tau*theta) * uValues(:,2) .* vValues(:,1);
elseif strcmp(dir, 'neg')
	T2 = sqrt(tau*theta) * uValues(:,3) .* vValues(:,2);
	T3 = sqrt(tau*theta) * uValues(:,1) .* vValues(:,1);
end
y = -T1 + T2 - T3;
end

function IMEXGK3(obj)
	electronConcentration = LegendrePoly(obj.X,  reshape(obj.coeffs,n+1,N), n);
	% get cell values
	CellValues = zeros(:,4);
	CellValues(:,2:3) = electronConcentration.getNodeValues';
	CellValues(:,1) = circshift(CellValues(:,3),1);
	CellValues(:,4) = circshift(CellValues(:,2),-1);
	
	% get qm1
	T1 = cellfun(@(f) gaussLegendre(@(x) 
end