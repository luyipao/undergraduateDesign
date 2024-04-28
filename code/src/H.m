function y = H(obj, Ecoeffs, ncoeffs)
rbbv = repmat(obj.basisBoundaryValues, obj.CellsNum, 1);
PPDP

EnFlux = 0.5 * (EValues(:,4) .* nValues(:,4) + EValues(:,3) .* nValues(:,3));
FEnv1 = mu * EnFlux .* rbbv(:,1);

EnFlux = 0.5 * (EValues(:,2).*nValues(:,2) + EValues(:,1).*nValues(:,1));
FEnv2 = mu * EnFlux .* rbbv(:,1);

for j = 1:N
	Aj = zeros(obj.degree+1);
	for i = 1:obj.degree+1
		for l = 1:obj.degree+1
		Aj(l,i) = PPDP(:,i,l) * Ecoeffs(:,j);
		end
	end
	A{j} = Aj;
end
A = - mu * blkdiag(A{:});

y = A * ncoeffs + FEnv1 - FEnv2;
end

% B only depend on basis functions values at each Cells.
function B = Hpos(obj)
B = repmat({obj.PnDPm}, 1, obj.CellsNum); % 将矩阵A重复n次并且形成一个数组
B = blkdiag(B{:});

D = obj.basisBoundaryValues(:,1) * obj.basisBoundaryValues(:,2)';
D = repmat({D}, obj.CellsNum, 1);
D = blkdiag(D{:});
D = circshift(D,-obj.degree-1,1);

E = obj.basisBoundaryValues(:,1) * obj.basisBoundaryValues(:,1)';
E = repmat({E}, obj.CellsNum, 1);
E = blkdiag(E);

B = B + D + E;
end

function B = Hneg(obj)
B = repmat({obj.PnDPm}, 1, obj.CellsNum); 
B = blkdiag(B{:});

D = obj.basisBoundaryValues(:,2) * obj.basisBoundaryValues(:,2)';
D = repmat({D}, obj.CellsNum, 1);
D = blkdiag(D{:});

E = obj.basisBoundaryValues(:,2) * obj.basisBoundaryValues(:,1)';
E = repmat({E}, obj.CellsNum, 1);
E = blkdiag(E);
E = circshift(E, obj.degree+1,1);

B = B + D + E;
end

function IMEXGK3(obj)
	electronConcentration = LegendrePoly(obj.X,  reshape(obj.coeffs,n+1,N), n);
	% get cell values
	CellValues = zeros(:,4);
	CellValues(:,2:3) = electronConcentration.getNodeValues';
	CellValues(:,1) = circshift(CellValues(:,3),1);
	CellValues(:,4) = circshift(CellValues(:,2),-1);
	
	[obj.Ecoeffs, obj.potentialCoeffs] = obj.getIMEXPotential;
	obj.H(obj.Eoceffs, obj.coeffs)
end


function [A, b] = IMEXPossionRelation(obj)
	% get relationship bettween field and potential CE = B * CP + E;
	I = eyes(obj.CellsNum);
	I = sparse(I);
	B = kron(obj.PnDPm, I);
	
	temp = obj.basisBoundaryValues(:,1) * obj.basisBoundaryValues(:,2)';
	I(1,1) = 0;
	D = kron(temp, I);
	D = circshift(D, -obj.degree-1, 1);
	
	E = obj.basisBoundaryValues(:,1) * obj.basisBoundaryValues(:,1);
	E = kron(temp, I);
	
	temp = zeros(obj.CellsNum,1);
	temp = sparse(temp);
	temp(end) = 1;
	b = -kron(1.5*obj.basisBoundaryValues(:,2), temp);
	
	A = B - D + E;
end

function [z, x] = getIMEXPotential(obj)
	% get relationship
	[A, b] = IMEXPossionRelation(obj);

	% get electric field z and electric potential  x
	I = eyes(obj.CellsNum);
	I = sparse(I);
	D = obj.basisBoundaryValues(:,1) * obj.basisBoundaryValues(:,1)';
	G = obj.basisBoundaryValues(:,1) * obj.basisBoundaryValues(:,2)';
	F = obj.basisBoundaryValues(:,2) * obj.basisBoundaryValues(:,2)';
	
	AE = kron(obj.PnDPm+D, I);
	BE = kron(G', diag(ones(obj.CellsNum-1,1), 1));
	CE = kron(G', diag(ones(obj.CellsNum-1, -1), -1);
	CE(1:obj.degree+1, 1:obj.degree+1) = F;
	
	AP = kron(G', sparse(diag(ones(obj.CellsNum-1,1), 1)));
	BP = kron(D+F, I);
	CP = kron(G,I);
	CP(1:obj.degree+1, 1:obj.degree+1) = zeros(obj.degree+1);
	dP = zeros(obj.CellsNum*(obj.degree+1),1);
	dP(1:obj.degree+1) = 1.5 * obj.basisBoundaryValues(:,2);
	
	x = ((AE-BE-CE)*A + AP-BP+CP ) / ((AE-BE-CE)*b - dP - 0.001546423010635 * (obj.coeffs - obj.initialCoeffs));
	z = A*x + b;
end