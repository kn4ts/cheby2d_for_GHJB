classdef ChebySeries2D
	properties
		M	% number of order
		Phi	% porinomial vector
		coef	% coefficient vector
		D1	% differential matrix with x1
		D2	% differential matrix with x2
		P0j	% product basis matrix
		Pi0	% product basis matrix
		Den	% denominator vector for integral calculation
		B	% product matrix (open-type)
		f	% function handle
	end
	methods
		function obj = ChebySeries2D(m,x)
			obj.M = m;
			obj.Phi = cheby2d_series(m);
			[obj.D1,  obj.D2 ] = genDifferentialMatrices_cheby2d_series(m);
			[obj.P0j, obj.Pi0] = genProductTensors_cheby2d_series(m);
			obj.Den = eval_innerproduct_cheby2d(m);
			if nargin >= 2
				if isa(x,'function_handle')
					obj = obj.approx(x);
				elseif isa(x,'double')
					obj.coef = x;
				end
			end
		end
		% ====
		% Calculate coef approximates function f(x1,x2)
		% ====
		function obj = approx(obj,f)
			N = obj.M^2;
			Num = [];
			w = @(x1,x2) 1./(sqrt(1-x1.^2).*sqrt(1-x2.^2)) ; % Weighting function for cheby
			for i=1:N
				integrated_fun = @(x1,x2) w(x1,x2) .* f(x1,x2) .* obj.Phi{i}(x1,x2) ;
				temp = integral2(integrated_fun, -1, 1, -1, 1);
				Num(i,1) = temp;
			end
			coef = Num./obj.Den;
			obj.coef = coef;
		end
		% ====
		% Calculate product of chebyshev porinomials
		% ====
		function obj3 = product(obj,obj2)
			%
			% Evaluate coefficient matrix for product of chebyshev 2d polynomials
			%  <a,Phi>*<b,Phi> = Te_x1' * C * Te_x2,  Te_x in R^(2M)
			%
			M = obj.M;
			%A = tensorize_2D_from_vector(a);
			A = tensorize_2D_from_vector(obj.coef);
			B = tensorize_2D_from_vector(obj2.coef);
			Be = [ B  zeros(M,M) ;
			       zeros(M,2*M)   ];
			temp = zeros(2*M,2*M);
			for i=1:M
				for j=1:M
					temp = temp + A(i,j)*obj.Pi0(:,:,i)*Be*obj.P0j(:,:,j);
				end
			end
			C = temp;
			%
			% Generate new chebyshev porinomial with coef of product result
			%
			obj3 = ChebySeries2D(M);
			obj3.coef = vectorize_from_2D_tensor(C(1:M,1:M));
			%obj3 = obj3.setCoefficient(vectorize_from_2D_tensor(C(1:M,1:M)));
			%obj3 = genOpenProductMatrix(obj3);
		end
		% ====
		% Evaluate open-type product matrix for 2D chebyshev polynomial series
		% B satisfies Phi * <b,Phi> = B * Phi
		% ====
		function B = productOpen(obj)
			%M = size(P0j,3);
			M = obj.M;
			M2 = M^2;
			B = [];
			a = ChebySeries2D(M);
			for i=1:M2
				a.coef = zeros(M2,1);
				a.coef(i,1) = 1;
				b = a.product(obj);
				%d_temp = vectorize_from_2D_tensor(Ctemp(1:M,1:M));
				B = [ B ; b.coef' ];
				%Popen = [ Popen ; d_temp' ];
			end
			obj.B = B;
		end
		function f = genFunc(obj)
			f = @(x1,x2) 0;
			for i=1:length(obj.Phi)
				f = @(x1,x2) f(x1,x2) + obj.coef(i) * obj.Phi{i}(x1,x2);
			end
		end
		function C = showTensolCoef(obj)
			C = tensorize_2D_from_vector(obj.coef);
		end
	end
end

% Generating 2D chebyshev polynomial series (M^2 functions)
%  from 1D chebyshev polynomial series T
% Phi(x1,x2) is R^(1 times M^2) vector defined as,
%  Phi = [ T_0(x1)T_0(x2)  T_0(x1)T_1(x2)  ... T_0(x1)T_M-1(x2) ...
%          T_1(x1)T_0(x2) ...              ... T_1(x1)T_M-1(x2) ...
%          ...
%          T_M-1(x1)T_0(x2) ...          ... T_M-1(x1)T_M-1(x2) ]
function Phi = cheby2d_series(M)
	T = cheby1d_series(M);
	%M = size(T,1);
	for i=1:M
		for j=1:M
			%Phi_temp = @(x) T{i}(x(1)).*T{j}(x(2));
			Phi_temp = @(x1,x2) T{i}(x1).*T{j}(x2);
			Phi{M*(i-1)+j} = Phi_temp;
		end
	end
	Phi = transpose(Phi);
end

% Generating 1D chebyshev polynomial series (M functions)
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
	if M == 1
		D1 = 0;
		D2 = 0;
	end
end

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

% Calculate inner prodct of 2D chebyshev function (series)
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

function V = tensorize_2D_from_vector(v)
	M2 = size(v,1);
	M  = sqrt(M2);
	V  = reshape(v,M,M)';
end

function v = vectorize_from_2D_tensor(V)
	M = size(V,1);
	v = reshape(V',M^2,1);
end

