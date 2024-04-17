function F = L(obj)
n = obj.degree;
N = obj.CellsNum;
parforX = obj.X;
parforCoeffs = reshape(obj.coeffs,n+1,N);


CellFuncs = cells{N}; 
CellValues = zeros(N,4);%[ll lr rl rr]

parforFunc2s = cells(N);
parforFluxs = zeros(N,1);

parfor j = 1:N
	% get electron concentration and cell values
	coeffsj = parforCoeffs(:,j);
	Cellj = obj.Cells(j);
	a = Cellj.a;
	b = Cellj.b;
	fProj = @(x) 0 * x;
	for i = 1:n+1
		fProj = @(x) fProj(x) + coeffsj(i) * Cellj.basisFunctions{i,1}(x);
	end
	CellFuncs{j} = fProj;

	fa = fProj(a);	
	fb = fProj(b);
	CellValues(j,2) = fa;
	CellValues(j,3) = fb;
	% get primitive electron concentration
	% input: obj.coeffs parforX
	% ouput: primitive electron concentration
	func2 = @(x) 0*x;
	for i = 1:n+1
		func2 = @(x) func2(x) + coeffsj(i) * sqrt((2*i-1)*h)/2 * obj.priLegendreFunctions{i}((2*x - 2*parforXc(j))/h);
	end
	parforFluxs(j) = func2(parforX(j+1)) - func2(parforX(j));
	parforFunc2s{j} = @(x) func2(x);
end

% get cell values
% 
CellValues(2:N-1, 1) = CellValues(1:N-2, 3);
CellValues(2:N-1, 4) = CellValues(3:N, 2);

% boundary 
CellValues(1, 1) = CellValues(N, 3);
CellValues(1, 4) = CellValues(2, 2);
CellValues(N, 1) = CellValues(N-1, 3);
CellValues(N, 4) = CellValues(1, 2);


% get primitive electron concentration
parforFluxs = [0; parforFluxs];
parforFluxs = cumsum(parforFluxs);
parfor j = 1:N
	parforFluxs(j) = -parforFunc2s{j}(parforX(j)) + parforFluxs(j);
end
for j = 1:N-1
	priElectronConcentration = @(x) priElectronConcentration(x) + (x>=parforX(j) & x < parforX(j+1)) .* (parforFunc2s{j}(x) + parforFluxs(j));
end
priElectronConcentration = @(x) priElectronConcentration(x) + (x>=parforX(N) & x <= parforX(N+1)) .* (parforFunc2s{N}(x) + parforFluxs(N));


% get electric field
T1 = gaussLegendre(@(x) priElectronConcentration(x), 0, 1);
T2 = gaussLegendre(@(x) obj.priDopingFunction(x), 0, 1); % set in generate function
electricField0 = 0.0015 * (T1 - T2);
electricField = @(x) -0.0015 * (priElectronConcentration(x) - obj.priDopingFunction(x) - priElectronConcentration(0) + priDopingFunction(0)) + electricField0 -  1.5;


% get auxiliary function 
% input: Cellj.func
auxCellValues = zeros(N,4); % [ll lr rl rr]
auxCellFuncs = cells(N,1);
auxCoeffs = zeros(n+1,N);
parfor j = 1:N
	Cellj = obj.Cells(j);
	temp = zeros(n+1,1);
	a = Cellj.a;
	b = Cellj.b;
	for i = 1:n+1
		tempf = @(x) Cellj.basisFunctions{i,2}(x) .* Cellj.func(x);
		temp(i) = - gaussLegendre(tempf,a,b);
	end
	
	Fj = CellValues(j,4) * Cellj.basisFunctionsBoundaryValues(:,2) - CellValues(j,2) * Cellj.basisFunctionsBoundaryValues(:,1) + temp;
	Fj =  0.139219332249189 * Fj;
	auxCoeffs(:,j) = Fj;
	
	auxq = @(x) 0*x;
	for i = 1:n+1
		auxq = @(x) auxq(x) + Fj(i) * Cellj.basisFunctions{i,1}(x);
	end
	auxCellFuncs{j} = auxq;

	auxCellValues(j,2) = auxq(parforX(j));
	auxCellValues(j,3) = auxq(parforX(j+1));
end

% 
auxCellValues(2:N-1, 1) = auxCellValues(1:N-2, 3);
auxCellValues(2:N-1, 4) = auxCellValues(3:N, 2);
% boundary
auxCellValues(1, 1) = auxCellValues(N, 3);
auxCellValues(1, 4) = auxCellValues(2, 2);
auxCellValues(N, 1) = auxCellValues(N-1, 3);
auxCellValues(N, 4) = auxCellValues(1, 2);

% L output F
F = zeros(n+1,N);
parfor j = 1:N
	Cellj = obj.Cells(j);
	a = Cellj.a;
	b = Cellj.b;
	Fj = F(:,j);
	parfor l = 1:n+1
		T1 = gaussLegendre(@(x) 0.75 * E(x) .* CellFuncs{j}(x) .* Cellj.basisFunctions{l,2}(x), a, b);
		T2 = gaussLegendre(@(x) 0.139219332249189 * auxCellFuncs{j}(x) .* Cellj.basisFunctions{l,2}(x), a, b);
		
		Eb = E(b);
		T3 = 0.75 * (max(Eb,0) * CellValues(j,4) + min(Eb,0) * CellValues(j,3));
		T3 = T3 + 0.139219332249189 * auxCellValues(j,3);
		T3 = T3 * Cellj.basisFunctionsBoundaryValues(l,2);
		
		Ea = E(a)
		T4 = 0.75 * (max(Ea,0) * CellValues(j,2) + min(Ea,0) * CellValues(j,1));
		T4 = T4 + 0.139219332249189 * auxCellValues(j,1);
		T4 = T4 * Cellj.basisFunctionsBoundaryValues(l,1);
		Fj(l) = -T1-T2+T3-T4;
	end
	F(:,j) = Fj;
end

end