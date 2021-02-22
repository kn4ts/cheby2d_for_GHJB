% Evaluate open-type product matrix for 2D chebyshev polynomial series
%  D satisfies Phi * <a,Phi> = D * Phi
function D = eval_OpenProductMatrix(b,P0j,Pi0)
	M = size(P0j,3);
	M2 = M^2;
	D = [];
	for i=1:M2
		a = zeros(M2,1);
		a(i,1) = 1;
		Ctemp = eval_CoefficientProductMatrix_cheby2d_series(a,b,P0j,Pi0);
		d_temp = vectorize_from_2D_tensor(Ctemp(1:M,1:M));
		D = [ D ; d_temp' ];
	end
end
