% Generating Differential Matrix for 2D chebyshev series
%  D1 satisfies ∂Phi/∂x1 = D1 * Phi, and
%  D2 satisfies ∂Phi/∂x2 = D2 * Phi
function [ D1, D2 ] = genDifferentialMatrices_cheby2d_series(M)
	R = zeros(M-1, M);
	R(1,1) = 1;
	for k=2:M-1
		r = 2*k ;
		%R(k,k) = 2*k;
		for j=k:-2:0
			if j>1
				R(k,j) = r;
			elseif j==1
				R(k,1) = 0.5*r;
			end
		end
	end
	R_fix = [ zeros(1,M) ; R ];
	D1 = kron(R_fix,eye(M)) ;
	D2 = kron(eye(M),R_fix) ;
end
