% Evaluate coefficient matrix for product of chebyshev 2d porinomials
%  <a,Phi>*<b,Phi> = Te_x1' * C * Te_x2,  Te_x in R^(2M)
function C = eval_CoefficientProductMatrix_cheby2d_series(a,b,P0j,Pi0)
	M = size(P0j,3);
	A = tensorize_2D_from_vector(a);
	B = tensorize_2D_from_vector(b);
	Be = [ B  zeros(M,M) ;
	       zeros(M,2*M)   ];
	temp = zeros(2*M,2*M);
	for i=1:M
		for j=1:M
			temp = temp + A(i,j)*Pi0(:,:,i)*Be*P0j(:,:,j);
		end
	end
	C = temp;
end
