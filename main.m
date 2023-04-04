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
M  = 7; % porinomial order
N  = 10; % iteration number

%x0 = [0.2;0];
x0 = zeros(2,1);

[ f1, f2, g1, g2, h1, h2 ] = genDynamics(sys) ;
[ f1L, f2L, g1L, g2L, A, B ] = genDynamicsLinear(sys,x0) ;
[K,S,e] = lqr(A,B,eye(2),1);
K

l  = @(x1,x2) 25*(x1.^2 +x2.^2);
u0 = @(x1,x2) -K(1)*x0(1) -K(2)*x0(2);

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


% for iteration
V_h = zeros(M^2,N+1);
V_h(:,1) = V_b.coef;

% iteration
for k=1:N

	dL_b = [ ChebySeries2D(M, Chv.D1'*V_b.coef) ;
		 ChebySeries2D(M, Chv.D2'*V_b.coef) ];
	% update input u
	GdL_b = [ G_b(1).product(dL_b(1));
		  G_b(2).product(dL_b(2)) ];
	u_b  = ChebySeries2D(M, -0.5*( GdL_b(1).coef +GdL_b(2).coef) );

	% evaluate ||u||^2
	uu_b = u_b.product(u_b);

	% evaluate product of matrices
	Gu_b = [ G_b(1).product(u_b);
		 G_b(2).product(u_b) ];

	% solve coefficient vector V
	E = 	Chv.D1 * F_b(1).productOpen + Chv.D1 * Gu_b(1).productOpen + ...
		Chv.D2 * F_b(2).productOpen + Chv.D2 * Gu_b(2).productOpen ;
	Vcoef = -pinv(E') * ( uu_b.coef + l_b.coef ) ;
	V_b = ChebySeries2D(M, Vcoef);

	% save coefficient
	V_h(:,k+1) = Vcoef;

end

%V_coef = (V);
U_coef = u_b.showTensolCoef

figure(1)
plot(V_h');

%u = genUfunc(Coef_u,Phi) ;
u = u_b.genFunc ;
uL = @(x1,x2) -K(1)*x1 -K(2)*x2;
V = V_b.genFunc ;

%dx = 0.01;
[X1,X2] = meshgrid(-1:0.1:1);
U  = u(X1,X2);
UL = uL(X1,X2);
L  = V(X1,X2);

figure(2)

subplot(221)
mesh(X1,X2,U)
title('input u (chebyshev porinomial)');
xlabel('z_1 (= 5x_1)');
ylabel('z_2 (= 5x_2)');

subplot(222)
mesh(X1,X2,UL)
title('input u_L (linear state feedback)');
xlabel('z_1 (= 5x_1)');
ylabel('z_2 (= 5x_2)');

subplot(223)
hold on
mesh(X1,X2,U)
mesh(X1,X2,UL)
title('comparison u and u_L');
xlabel('z_1 (= 5x_1)');
ylabel('z_2 (= 5x_2)');
grid on
hold off

subplot(224)
mesh(X1,X2,L)
title('function V_b');
xlabel('z_1 (= 5x_1)');
ylabel('z_2 (= 5x_2)');


% ======================================
x0 = [2.5;0];
Ts = 0.01;
Tend = 20;
Ns = Tend/Ts;
TimeSeries = 0:Ts:Tend-Ts ;

for Flag=0:1

	x = x0;
	x_h = zeros(2,Ns);
	u_h = zeros(1,Ns);

	x_h(:,1) = x;

	for i=1:Ns-1
	
		switch Flag
		case 0
			uv = uL(0.2*x(1),0.2*x(2));
		case 1
			uv = u(0.2*x(1),0.2*x(2));
		end
	
		xd =  [ f1(0.2*x(1),0.2*x(2))+g1(0.2*x(1),0.2*x(2))*uv ;
			f2(0.2*x(1),0.2*x(2))+g2(0.2*x(1),0.2*x(2))*uv ];
		x  = x + Ts * xd ;

		x_h(:,i+1) = x;
		u_h(:,i+1) = uv;
	end

	figure(3)

	subplot(2,1,1)
	hold on
	switch Flag
		case 0
			plot(TimeSeries,x_h,'--');
		case 1
			plot(TimeSeries,x_h);
	end
	title('state variables');
	grid on

	subplot(2,1,2)
	hold on
	switch Flag
		case 0
			plot(TimeSeries,u_h,'--');
		case 1
			plot(TimeSeries,u_h);
	end
	title('input signal');
	%plot(TimeSeries,u_h);
	grid on

end


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

	%f1 = @(x1,x2) x2 ;
	%f2 = @(x1,x2) a1*sin(x1) + a2*x2 ;
	%g1 = @(x1,x2) 0 ;
	%g2 = @(x1,x2) b1 ;
	%h1 = @(x1,x2) sys.l * sin(x1) ;
	%h2 = @(x1,x2) sys.l * (1-cos(x1)) ;

	f1 = @(x1,x2) 5*x2 ;
	f2 = @(x1,x2) a1*sin(5*x1) + a2*5*x2 ;
	g1 = @(x1,x2) 0 ;
	g2 = @(x1,x2) b1 ;
	h1 = @(x1,x2) sys.l * sin(5*x1) ;
	h2 = @(x1,x2) sys.l * (1-cos(5*x1)) ;
end

function [ f1, f2, g1, g2, A, B ] = genDynamicsLinear(sys, x0)
	%z0 = 5*x0 ;
	A = [ 0  1;
	      %-sys.g*cos(x0(1,:))/sys.l  -sys.mu ];
	      -sys.g*cos(x0(1,:))/sys.l  -sys.mu/(sys.m*sys.l^2) ];
	B = [ 0 ;
	      1/(sys.m*sys.l^2) ];
	f1 = @(x1,x2) A(1,1)*x1 +A(1,2)*x2 ;
	f2 = @(x1,x2) A(2,1)*x1 +A(2,2)*x2 ;
	g1 = @(x1,x2) B(1,1) ;
	g2 = @(x1,x2) B(2,1) ;
end
