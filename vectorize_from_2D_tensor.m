function v = vectorize_from_2D_tensor(V)
	M = size(V,1);
	v = reshape(V',M^2,1);
end
