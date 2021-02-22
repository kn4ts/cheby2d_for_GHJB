% Caliculate inner prodct of 2D chebyshev function (series)
function Den = eval_innerproduct_cheby2d(M)
	%N = size(Phi,2);
	N = M^2;
	for i=1:N
		if i==1
			temp = pi^2;
		elseif i<=M
			temp = 0.5*pi^2;
		elseif rem(i,M)==1
			temp = 0.5*pi^2;
		else
			temp = 0.25*pi^2;
		end
		Den(i,1) = temp;
	end
end
