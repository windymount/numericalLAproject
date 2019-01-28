function [P,U,V] = uzawa(P,U,V,r,itertimes)
	N = length(P);
	res = zeros(N-1,2*N);
	for j = 1:N
		for i = 1:N-1
			res(i,j) = (P(i+1,j)-P(i,j))*N;
			res(i,j+N) = (P(j,i+1)-P(j,i))*N;
		end
    end
    U = [zeros(1,N); U; zeros(1,N)];
	V = [zeros(1,N); V; zeros(1,N)];
	tr = r - res;
	for m = 1:itertimes
		for i = 1:N
			for j = 2:N
				if i== 1
					U(j,i) = (N^2*(U(j,i+1)+U(j+1,i)+U(j-1,i))+tr(j-1,i))/(5 * N^2);
					V(j,i) = (N^2*(V(j,i+1)+V(j+1,i)+V(j-1,i))+tr(j-1,N+i))/(5 * N^2);
				elseif i== N
					U(j,i) = (N^2*(U(j,i-1)+U(j+1,i)+U(j-1,i))+tr(j-1,i))/(5 * N^2);
					V(j,i) = (N^2*(V(j,i-1)+V(j+1,i)+V(j-1,i))+tr(j-1,N+i))/(5 * N^2);
				else
					U(j,i) = (N^2*(U(j,i+1)+U(j,i-1)+U(j+1,i)+U(j-1,i))+tr(j-1,i))/(4 * N^2);
					V(j,i) = (N^2*(V(j,i+1)+V(j,i-1)+V(j+1,i)+V(j-1,i))+tr(j-1,N+i))/(4 * N^2);
				end
			end
		end
    end
	for i = 1:N
		for j = 1:N
			P(i,j)= P(i,j)-0.03*((U(i+1,j) - U(i,j)) +(V(j+1,i)-V(j,i)))*N;
		end
	end
	U = U(2:N,:);
	V = V(2:N,:);