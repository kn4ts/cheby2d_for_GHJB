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
)

x0 = zeros(2,1);

[ f1, f2, g1, g2, h1, h2 ] = genDynamics(sys) ;
%[ f1, f2, g1, g2 ] = genDynamicsLinear(sys,x0) ;
l = @(x1,x2) x1.^2 +x2.^2;
%P = [ 0.05  0 ;
%      0   0.5 ];

u0  = @(x1,x2) 0;

% Number of terms of Chebyshev Porinomial
M = 5;
N  = 10;

T   = cheby1d_series(M);
Phi = cheby2d_series(T);

Den = eval_innerproduct_cheby2d(M);

% Prepare operation matrices for chebyshev porinomial series
[ D1, D2 ] = genDifferentialMatrices_cheby2d_series(M);
[ P0j, Pi0 ] = genProductTensors_cheby2d_series(M);

% Fit functions
Numf1 = eval_integral2_fun_w_cheby2d(f1,Phi);
Numf2 = eval_integral2_fun_w_cheby2d(f2,Phi);

Numg1 = eval_integral2_fun_w_cheby2d(g1,Phi);
Numg2 = eval_integral2_fun_w_cheby2d(g2,Phi);

Numl  = eval_integral2_fun_w_cheby2d(l,Phi);

Num_u0  = eval_integral2_fun_w_cheby2d(u0,Phi);

Coef_f1 = Numf1./Den;
Coef_f2 = Numf2./Den;

Coef_g1 = Numg1./Den;
Coef_g2 = Numg2./Den;

Coef_u0  = Num_u0./Den;

Coef_l = Numl./Den;

% Calculate coefficient matrices
C_G1u = eval_CoefficientProductMatrix_cheby2d_series(Coef_g1, Coef_u0, P0j, Pi0);
C_G2u = eval_CoefficientProductMatrix_cheby2d_series(Coef_g2, Coef_u0, P0j, Pi0);

Coef_G1u = vectorize_from_2D_tensor(C_G1u(1:M,1:M));
Coef_G2u = vectorize_from_2D_tensor(C_G2u(1:M,1:M));

C_uu  = eval_CoefficientProductMatrix_cheby2d_series(Coef_u0, Coef_u0, P0j, Pi0);
Coef_un  = vectorize_from_2D_tensor(C_uu(1:M,1:M));

F1  = eval_OpenProductMatrix(Coef_f1, P0j, Pi0);
F2  = eval_OpenProductMatrix(Coef_f2, P0j, Pi0);

G1u = eval_OpenProductMatrix(Coef_G1u, P0j, Pi0);
G2u = eval_OpenProductMatrix(Coef_G2u, P0j, Pi0);

% Solve linear first order equation
E = D1*F1 +D1*G1u +D2*F2 +D2*G2u;
V = -pinv(E')*(Coef_un+Coef_l);

% for iteration
V_h = zeros(M^2,N+1);
V_h(:,1) = V;

% iteration
for k=1:N

	% update input u
	C_g1dL = eval_CoefficientProductMatrix_cheby2d_series(Coef_g1, D1'*V, P0j, Pi0);
	C_g2dL = eval_CoefficientProductMatrix_cheby2d_series(Coef_g2, D2'*V, P0j, Pi0);
	C_u    = -0.5* (C_g1dL +C_g2dL) ;
	Coef_u = vectorize_from_2D_tensor(C_u(1:M,1:M));

	% evaluate ||u||^2
	C_un    = eval_CoefficientProductMatrix_cheby2d_series(Coef_u, Coef_u, P0j, Pi0);
	Coef_un = vectorize_from_2D_tensor(C_un(1:M,1:M));

	% evaluate product of matrices
	C_G1u = eval_CoefficientProductMatrix_cheby2d_series(Coef_g1, Coef_u, P0j, Pi0);
	C_G2u = eval_CoefficientProductMatrix_cheby2d_series(Coef_g2, Coef_u, P0j, Pi0);
	Coef_G1u = vectorize_from_2D_tensor(C_G1u(1:M,1:M));
	Coef_G2u = vectorize_from_2D_tensor(C_G2u(1:M,1:M));

	G1u = eval_OpenProductMatrix(Coef_G1u, P0j, Pi0);
	G2u = eval_OpenProductMatrix(Coef_G2u, P0j, Pi0);

	% solve coefficient vector V
	E = D1*F1 +D1*G1u +D2*F2 +D2*G2u;
	V = -pinv(E')*(Coef_un+Coef_l);

	% save coefficient
	V_h(:,k+1) = V;

end

V_coef = tensorize_2D_from_vector(V)
U_coef = tensorize_2D_from_vector(Coef_u)
plot(V_h');

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

function [ f1, f2, g1, g2 ] = genDynamicsLinear(sys, x0)
	A = [ 0  1;
	      -sys.g*cos(x0(1,:))/sys.l  -sys.mu ];
	B = [ 0 ;
	      1/(sys.m*sys.l^2) ];
	f1 = @(x1,x2) A(1,1)*x1 +A(1,2)*x2 ;
	f2 = @(x1,x2) A(2,1)*x1 +A(2,2)*x2 ;
	g1 = @(x1,x2) B(1,1) ;
	g2 = @(x1,x2) B(2,1) ;
end
