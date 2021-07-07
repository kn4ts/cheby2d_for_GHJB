clear all
close all
clc

profile on

% System in paper [1]
% [1] Approximate Solutions to Hamilton-Jacobi Equations Based on Chebyshev Polynomials,
%     URL: https://www.jstage.jst.go.jp/article/sicetr1965/44/2/44_133/_article
sys = struct( ...
	'm', 1, ... % kg
	'l', 1, ... % m
	'mu', 1, ...% Ns
	'g', 9.8 ...% m/s^2
);

% Number of terms of Chebyshev Porinomial
M = 10;
N  = 10;

x0 = zeros(2,1);

[ f1, f2, g1, g2, h1, h2 ] = genDynamics(sys) ;
[ f1L, f2L, g1L, g2L, A, B ] = genDynamicsLinear(sys,x0) ;
l  = @(x1,x2) x1.^2 +x2.^2;
u0 = @(x1,x2) 0;

Chv = ChebySeries2D(M);
F_b = [ ChebySeries2D(M, f1);
	ChebySeries2D(M, f2) ];
G_b = [ ChebySeries2D(M, g1);
	ChebySeries2D(M, g2) ];
l_b =   ChebySeries2D(M, l );
u_b =   ChebySeries2D(M, u0);

Gu_b = [ G_b(1).product(u_b);
	 G_b(2).product(u_b) ];
uu_b = u_b.product(u_b);

E = 	Chv.D1 * F_b(1).productOpen + Chv.D1 * Gu_b(1).productOpen + ...
	Chv.D2 * F_b(2).productOpen + Chv.D2 * Gu_b(2).productOpen ;
Vcoef = -pinv(E') * ( uu_b.coef + l_b.coef ) ;
V_b = ChebySeries2D(M, Vcoef);
%V_b = ChebySeries2D(M);
%V_b.coef = -pinv(E') * ( uu_b.coef + l_b.coef ) ;
%f1f2_b = f1_b.product(f2_b)
%      0   0.5 ];
[K,S,e] = lqr(A,B,eye(2),1);
K

% for iteration
V_h = zeros(M^2,N+1);
V_h(:,1) = V_b.coef;

% iteration
for k=1:N

	%dL_b = [ ChebySeries2D(M) ;
	%	 ChebySeries2D(M) ];
	dL_b = [ ChebySeries2D(M, Chv.D1'*V_b.coef) ;
		 ChebySeries2D(M, Chv.D2'*V_b.coef) ];
	% update input u
	GdL_b = [ G_b(1).product(dL_b(1));
		  G_b(2).product(dL_b(2)) ];
	u_b  = ChebySeries2D(M, -0.5 * ( GdL_b(1).coef +GdL_b(2).coef) );
	%C_u    = -0.5* (C_g1dL +C_g2dL) ;
	%Coef_u = vectorize_from_2D_tensor(C_u(1:M,1:M));

	% evaluate ||u||^2
	uu_b = u_b.product(u_b);
	%C_un    = eval_CoefficientProductMatrix_cheby2d_series(Coef_u, Coef_u, P0j, Pi0);
	%Coef_un = vectorize_from_2D_tensor(C_un(1:M,1:M));

	% evaluate product of matrices
	Gu_b = [ G_b(1).product(u_b);
		 G_b(2).product(u_b) ];

	% solve coefficient vector V
	%E = D1*F1 +D1*G1u +D2*F2 +D2*G2u;
	%V = -pinv(E')*(Coef_un+Coef_l);
	E = 	Chv.D1 * F_b(1).productOpen + Chv.D1 * Gu_b(1).productOpen + ...
		Chv.D2 * F_b(2).productOpen + Chv.D2 * Gu_b(2).productOpen ;
	Vcoef = -pinv(E') * ( uu_b.coef + l_b.coef ) ;
	V_b = ChebySeries2D(M, Vcoef);
	% save coefficient
	V_h(:,k+1) = Vcoef;

end

%V_coef = (V);
U_coef = u_b.showTensolCoef

figure(2)
plot(V_h');

%u = genUfunc(Coef_u,Phi) ;
u = u_b.genFunc ;
uL = @(x1,x2) -K(1)*x1 -K(2)*x2;

%dx = 0.01;
[X1,X2] = meshgrid(-1:0.1:1);
U  = u(X1,X2);
UL = uL(X1,X2);

figure(1)
subplot(121)
mesh(X1,X2,U)
subplot(122)
mesh(X1,X2,UL)
%for i=1:200
%	for j=1:200
%		x_1 = dx * i - 1;
%		x_2 = dx * j - 1;
%		u_h(i,j) = x_1
%	end
%end

profile viewer


% System in paper[1]
% [1] Approximate Solutions to Hamilton-Jacobi Equations Based on Chebyshev Polynomials,
%     URL: https://www.jstage.jst.go.jp/article/sicetr1965/44/2/44_133/_article
function [ f1, f2, g1, g2, h1, h2 ] = genDynamics(sys)
	a1 = - sys.g/sys.l ;
	a2 = - sys.mu/(sys.m*sys.l^2) ;

	b1 = 1/(sys.m*sys.l^2) ;

	f1 = @(x1,x2) x2 ;
	f2 = @(x1,x2) a1*sin(x1) + a2*x2 ;
	g1 = @(x1,x2) 0 ;
	g2 = @(x1,x2) b1 ;
	h1 = @(x1,x2) sys.l * sin(x1) ;
	h2 = @(x1,x2) sys.l * (1-cos(x1)) ;
end
function [ f1, f2, g1, g2, h1, h2 ] = genDynamics_(sys)
	a1 = - sys.g/sys.l ;
	a2 = - sys.mu/(sys.m*sys.l^2) ;

	b1 = 1/(sys.m*sys.l^2) ;

	f1 = @(x) x(2) ;
	f2 = @(x) a1*sin(x(1)) + a2*x(2) ;
	g1 = @(x) 0 ;
	g2 = @(x) b1 ;
	h1 = @(x) sys.l * sin(x(1)) ;
	h2 = @(x) sys.l * (1-cos(x(1))) ;
end

function [ f1, f2, g1, g2, A, B ] = genDynamicsLinear(sys, x0)
	A = [ 0  1;
	      -sys.g*cos(x0(1,:))/sys.l  -sys.mu ];
	B = [ 0 ;
	      1/(sys.m*sys.l^2) ];
	f1 = @(x1,x2) A(1,1)*x1 +A(1,2)*x2 ;
	f2 = @(x1,x2) A(2,1)*x1 +A(2,2)*x2 ;
	g1 = @(x1,x2) B(1,1) ;
	g2 = @(x1,x2) B(2,1) ;
end
