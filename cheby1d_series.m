% Generating 1D chebyshev porinomial series (M functions)
% T = [ T_0(x)  T_1(x)  ... T_M(x) ]^T
function T = cheby1d_series(M)
	T = {};
	for i=0:M-1
		if i==0
			T_temp = @(x) 1;
		elseif i==1
			T_temp = @(x) x;
		else
			T_temp = @(x) 2.*x.*T{i}(x) -T{i-1}(x);
		end
		T{i+1} = T_temp ;
	end
	T = transpose(T);
end
