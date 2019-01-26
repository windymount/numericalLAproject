function [P,U,V] = DGS(P,U,V,r)
	N = length(P);
	%add zeros boundary
	U = [zeros(1,N); U; zeros(1,N) ];
	V = [zeros(1,N); V; zeros(1,N) ];
	% step 1 G-S iteration
	for i = 1:N
		for j = 2:N
			if i== 1 
				U(j,i) = (N^2*(U(j,i+1)+U(j+1,i)+U(j-1,i))-N *(P(j,i)-P(j-1,i))+r(j-1,i))/(3 * N^2);
				V(j,i) = (N^2*(V(j,i+1)+V(j+1,i)+V(j-1,i))-N *(P(i,j)-P(i,j-1))+r(j-1,N+i))/(3 * N^2);
			elseif i== N
				U(j,i) = (N^2*(U(j,i-1)+U(j+1,i)+U(j-1,i))-N *(P(j,i)-P(j-1,i))+r(j-1,i))/(3 * N^2);
				V(j,i) = (N^2*(V(j,i-1)+V(j+1,i)+V(j-1,i))-N *(P(i,j)-P(i,j-1))+r(j-1,N+i))/(3 * N^2);
			else
				U(j,i) = (N^2*(U(j,i+1)+U(j,i-1)+U(j+1,i)+U(j-1,i))-N *(P(j,i)-P(j-1,i))+r(j-1,i))/(4 * N^2);
				V(j,i) = (N^2*(V(j,i+1)+V(j,i-1)+V(j+1,i)+V(j-1,i))-N *(P(i,j)-P(i,j-1))+r(j-1,N+i))/(4 * N^2);
			end
		end
	end
	% step 2 calculate divergence
	div = zeros(N);
	for i = 1:N
		for j = 1:N
			div(i,j) = -((U(i+1,j) - U(i,j)) +(V(j+1,i)-V(j,i)))*N;
		end
	end
	% step 3 update P,V,U of inner points
	U(2:N-1,2:N-1) = U(2:N-1,2:N-1) - div(2:N-1,2:N-1)/(4*N);
	U(3:N,2:N-1)   = U(3:N,2:N-1)   + div(2:N-1,2:N-1)/(4*N);
	V(2:N-1,2:N-1) = V(2:N-1,2:N-1) - div(2:N-1,2:N-1)'/(4*N);
	V(3:N,2:N-1)   = V(3:N,2:N-1)   + div(2:N-1,2:N-1)'/(4*N);
	P(2:N-1,2:N-1) = P(2:N-1,2:N-1) + div(2:N-1,2:N-1);
	P(2:N-1,1:N-2) = P(2:N-1,1:N-2) - div(2:N-1,2:N-1)/4;
	P(2:N-1,3:N)   = P(2:N-1,3:N)   - div(2:N-1,2:N-1)/4;
	P(1:N-2,2:N-1) = P(1:N-2,2:N-1) - div(2:N-1,2:N-1)/4;
	P(3:N,2:N-1)   = P(3:N,2:N-1)   - div(2:N-1,2:N-1)/4;
	% step 4 update P,V,U on bonudary
	U(2:N-1,[1 N]) = U(2:N-1,[1 N]) - div(2:N-1,[1 N] )/(3*N);
	U(3:N,[1 N])   = U(3:N,[1 N]) + div(3:N,[1 N] )/(3*N);
	V(2:N-1,[1 N]) = V(2:N-1,[1 N]) - div([1 N],2:N-1)'/(3*N);
	V(3:N,[1 N])   = V(3:N,[1 N]) + div([1 N],3:N)'/(3*N);
	U(2,2:N-1) = U(2,2:N-1) + div(1,2:N-1)/(3*N);
	U(N,2:N-1) = U(N,2:N-1) - div(N,2:N-1)/(3*N);
	V(2,2:N-1) = V(2,2:N-1) + div(2:N-1,1)'/(3*N);
	V(N,2:N-1) = V(N,2:N-1) - div(2:N-1,N)'/(3*N);
	ker = [0 -1/3 0;-1/3 4/3 -1/3;0 -1/3 0];
	temp = div;
	temp([1 N],[1 N]) = zeros(2);
	temp(2:N-1,2:N-1) = zeros(N-2);
	temp = conv2(temp,ker);
	temp = temp(2:N+1,2:N+1);
	P = P + temp;
	U(2,[1 N]) = U(2,[1 N]) + div(1,[1 N])/(2*N);
	U(N,[1 N]) = U(N,[1 N]) - div(N,[1 N])/(2*N);
	V(2,[1 N]) = V(2,[1 N]) + div([1 N],1)'/(2*N);
	V(N,[1 N]) = V(N,[1 N]) - div([1 N],N)'/(2*N);
	temp = zeros(N);
	temp([1 N],[1 N]) = div([1 N],[1 N]);
	ker = [0 -0.5 0;-0.5 2 -0.5;0 -0.5 0];
	temp = conv2(temp,ker);
	temp = temp(2:N+1,2:N+1);
	P = P + temp;
	U = U(2:N,:);
	V = V(2:N,:);


