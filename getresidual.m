function r = getresidual(P,U,V)
	N = length(P);
	U = [zeros(1,N); U; zeros(1,N) ];
	V = [zeros(1,N); V; zeros(1,N) ];
	res = zeros(N-1,2*N);
	for j = 1:N
		for i = 1:N-1
			if j == 1
				res(i,j) = -(-3*U(i+1,j)+U(i,j)+U(i+2,j)+U(i+1,j+1))*N^2+(P(i+1,j)-P(i,j))*N;
				res(i,j+N) = -(-3*V(i+1,j)+V(i,j)+V(i+2,j)+V(i+1,j+1))*N^2+(P(j,i+1)-P(j,i))*N;
			elseif j == N
				res(i,j) = -(-3*U(i+1,j)+U(i,j)+U(i+2,j)+U(i+1,j-1))*N^2+(P(i+1,j)-P(i,j))*N;
				res(i,j+N) = -(-3*V(i+1,j)+V(i,j)+V(i+2,j)+V(i+1,j-1))*N^2+(P(j,i+1)-P(j,i))*N;
			else
				res(i,j) = -(-4*U(i+1,j)+U(i,j)+U(i+2,j)+U(i+1,j+1)+U(i+1,j-1))*N^2+(P(i+1,j)-P(i,j))*N;
				res(i,j+N) = -(-4*V(i+1,j)+V(i,j)+V(i+2,j)+V(i+1,j+1)+V(i+1,j-1))*N^2+(P(j,i+1)-P(j,i))*N;
			end
		end
	end
	r = res;