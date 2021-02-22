function V = tensorize_2D_from_vector(v)
	M2 = size(v,1);
	M  = sqrt(M2);
	V  = reshape(v,M,M)';
end
