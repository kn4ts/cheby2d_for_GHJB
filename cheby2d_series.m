% Generating 2D chebyshev porinomial series (M^2 functions)
%  from 1D chebyshev porinomial series T
% Phi(x1,x2) is R^(1 times M^2) vector defined as,
%  Phi = [ T_0(x1)T_0(x2)  T_0(x1)T_1(x2)  ... T_0(x1)T_M-1(x2) ...
%          T_1(x1)T_0(x2) ...              ... T_1(x1)T_M-1(x2) ...
%          ...
%          T_M-1(x1)T_0(x2) ...          ... T_M-1(x1)T_M-1(x2) ]
function Phi = cheby2d_series(T)
	M = size(T,1);
	for i=1:M
		for j=1:M
			Phi_temp = @(x1,x2) T{i}(x1).*T{j}(x2);
			Phi{M*(i-1)+j} = Phi_temp;
		end
	end
	Phi = transpose(Phi);
end
