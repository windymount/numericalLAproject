function [P,U,V] = multigrid(v1,v2,iterfunc,N)
	% N:number of grids in discretization
	% iterfunc:磨光子函数
	TOL = 1e-8;
	f = @(x,y)(-4*pi^2*(2*cos(2*pi*x)-1)*sin(2*pi * y)+x^2);
	g = @(x,y)(4*pi^2*(2*cos(2*pi*y)-1)*sin(2*pi * x));
	P = zeros(N);
	U = zeros(N-1,N);
	V = zeros(N-1,N);
	Pcell = cell(gridnum,1);
	Ucell = cell(gridnum,1);
	Vcell = cell(gridnum,1);
	rcell = cell(gridnum,1);
	% initialize residual
	% r:residual of Stokes function.Stored by a 2N*N-1 matrix
	r = zeros(N-1,2*N);
	for i = 1:N-1
		for j = 1:N
			r(i,j) = f(i/N,(j-0.5)*N);
		end
	end
	for i = 1:N-1
		for j = N+1:2*N
			r(i,j) = g((j-0.5)/N,i/N);
		end
	end
	rcell{1} = r;
	% V-cycle multigrid
	while True
		for i = 2:gridnum
			r = confine(r);
			rcell(i) = r;
			P = zeros(N/(2^i));
			U = zeros(N/(2^i)-1,N/(2^i));
			V = zeros(N/(2^i)-1,N/(2^i));
			for j = 1:v1
				[P,U,V] = iterfunc(P,U,V,r);
			end
			r = r - getresidual(P,U,V);
			Pcell(i) = P;
			Ucell(i) = U;
			Vcell(i) = V;
		end
		for i = gridnum-1:-1:1
			r = rcell{i};
			P = raiseP(Pcell{i+1})+Pcell{i};
			U = raiseU(Ucell{i+1})+Ucell{i};
			V = raiseU(Vcell{i+1})+Vcell{i};
			for j = 1:v2
				[P,U,V] = iterfunc(P,U,V,r);
			Pcell{i} = P;
			Ucell{i} = U;
			Vcell{i} = V;
			end
		end
		r = rcell{1} - getresidual(P,U,V);
		if norm(r) < TOL
			break
		end
	end
end
function r = confine(r)
	N = length(r)/2;
	ker = [0.125 0.25 0.125;0.125 0.25 0.125;0 0 0];
	r = conv2(r,ker);
	r = r(2:2:2*N,2:2:N-1);
end
function P = raiseP(P)
	N = length(P);
	ker = [0.25 0.25 0;0.25 0.25 0;0 0 0];
	P = conv2(P,ker);
	P = P(2:2:N,2:2:N);
end
function U = raiseU(U)
	N = length(U) + 1;
	ker = [0.125 0.25 0.125;0.125 0.25 0.125;0 0 0];
	U = conv2(U,ker);
	U = U(2:2:N,2:2:N-1);
end
