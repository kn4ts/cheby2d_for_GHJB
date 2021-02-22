% Generating Coefficient Matrix for producting 2D chebyshev function series
%  <a,Phi><b,Phi> = Te_x1' * C * Te_x2,  Te in R^(2M)
%  C = Σ_(i=1)^(M-1) Σ_(j=1)^(M-1) aij * Pi0 * Be * P0j
%  Pi0 = P0i' for i = 0 ~ M
function [ P0j, Pi0 ] = genProductTensors_cheby2d_series(M)
	P0j(:,:,1) = eye(2*M,2*M);
	Pi0(:,:,1) = eye(2*M,2*M);
	for i=1:M-1
		Temp = [];
		for k=1:M
			temp = zeros(1,2*M);
			if k==1
				temp(:,i+1) = 1;
			else
				if i-k+1>0
					temp(:,i-k+2) = 0.5;
				else
					temp(:,abs(i-k)) = 0.5;
				end
				temp(:,i+k) = 0.5;
			end
		Temp = [ Temp ; temp ];
		end
		P0j(:,:,i+1) = [ Temp ;
		               zeros(M,2*M) ];
		Pi0(:,:,i+1) = [ Temp ;
		               zeros(M,2*M) ]';
	end
end
