% Caliculate integral of f(x1,x2) with weight and 2D chebyshev function (series)
function Num = eval_integral2_fun_w_cheby2d(f,Phi)
	N = size(Phi,1);
	Num = [];
	w = @(x1,x2) 1./(sqrt(1-x1.^2).*sqrt(1-x2.^2)) ; % Weighting function
	for i=1:N
		integrated_fun = @(x1,x2) w(x1,x2) .* f(x1,x2) .* Phi{i}(x1,x2) ;
		temp = integral2(integrated_fun, -1, 1, -1, 1);
		Num(i,1) = temp;
	end
end
