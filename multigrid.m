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
	% r:residual of Stokes function
	r = zeros(2*N,N-1);
	for i = 1:N
		for j = 1:N-1
			r(i,j) = f((j-1)/N,(i-0.5)/N);
		end
	end
	for i = N+1:2*N
		for j = 1:N-1
			r(i,j) = g((i-0.5)/N,(j-1)/N);
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
