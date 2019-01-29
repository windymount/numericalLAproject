function [P,U,V] = DGS(P,U,V,r)
	N = length(P);
	%add zeros boundary
	U = [zeros(1,N); U; zeros(1,N)];
	V = [zeros(1,N); V; zeros(1,N)];
	% step 1 G-S iteration
	for i = 1:N
		for j = 2:N
			if i== 1 
				U(j,i) = (N^2*(U(j,i+1)+U(j+1,i)+U(j-1,i))-N *(P(j,i)-P(j-1,i))+r(j-1,i))/(5 * N^2);
				V(j,i) = (N^2*(V(j,i+1)+V(j+1,i)+V(j-1,i))-N *(P(i,j)-P(i,j-1))+r(j-1,N+i))/(5 * N^2);
			elseif i== N
				U(j,i) = (N^2*(U(j,i-1)+U(j+1,i)+U(j-1,i))-N *(P(j,i)-P(j-1,i))+r(j-1,i))/(5 * N^2);
				V(j,i) = (N^2*(V(j,i-1)+V(j+1,i)+V(j-1,i))-N *(P(i,j)-P(i,j-1))+r(j-1,N+i))/(5 * N^2);
			else
				U(j,i) = (N^2*(U(j,i+1)+U(j,i-1)+U(j+1,i)+U(j-1,i))-N *(P(j,i)-P(j-1,i))+r(j-1,i))/(4 * N^2);
				V(j,i) = (N^2*(V(j,i+1)+V(j,i-1)+V(j+1,i)+V(j-1,i))-N *(P(i,j)-P(i,j-1))+r(j-1,N+i))/(4 * N^2);
			end
		end
    end
	% step 2 update P,V,U of inner points
	for i = 2:N-1
		for j = 2:N-1
			div = -((U(i+1,j) - U(i,j)) +(V(j+1,i)-V(j,i)))*N;
			delta = div / 4 / N;
			U(i+1,j) = U(i+1,j) + delta;
			U(i,j) = U(i,j) - delta;
			V(j+1,i) = V(j+1,i) + delta;
			V(j,i) = V(j,i) - delta;
			P(i,j) = P(i,j) + div;
			P(i+1,j) = P(i+1,j) - div/4;
			P(i-1,j) = P(i-1,j) - div/4;
			P(i,j+1) = P(i,j+1) - div/4;
			P(i,j-1) = P(i,j-1) - div/4;
		end
	end
	% step 3 update P,V,U on bonudary
	for i = [1 N]
		for j = 2:N-1
			div = -((U(i+1,j) - U(i,j)) +(V(j+1,i)-V(j,i)))*N;
			delta = div / 3 / N;
			V(j+1,i) = V(j+1,i) + delta;
			V(j,i) = V(j,i) - delta;
			P(i,j) = P(i,j) + 4*div/3;
			P(i,j+1) = P(i,j+1) - div/3;
			P(i,j-1) = P(i,j-1) - div/3;
			if i == 1
				U(i+1,j) = U(i+1,j) + delta;
				P(i+1,j) = P(i+1,j) - div/3;
			else
				U(i,j) = U(i,j) - delta;
				P(i-1,j) = P(i-1,j) - div/3;
			end
		end
	end
	for j = [1 N]
		for i = 2:N-1
			div = -((U(i+1,j) - U(i,j)) +(V(j+1,i)-V(j,i)))*N;
			delta = div / 3 / N;
			U(i+1,j) = U(i+1,j) + delta;
			U(i,j) = U(i,j) - delta;
			P(i,j) = P(i,j) + 4*div/3;
			P(i+1,j) = P(i+1,j) - div/3;
			P(i-1,j) = P(i-1,j) - div/3;
			if j == 1
				V(j+1,i) = V(j+1,i) + delta;
				P(i,j+1) = P(i,j+1) - div/3;
			else
				V(j,i) = V(j,i) - delta;
				P(i,j-1) = P(i,j-1) - div/3;
			end
		end
	end
	for i = [1 N]
		for j = [1 N]
			div = -((U(i+1,j) - U(i,j)) +(V(j+1,i)-V(j,i)))*N;
			delta = div / 2 / N;		
			if i == 1 && j == 1
				U(i+1,j) = U(i+1,j) + delta;
				V(j+1,i) = V(j+1,i) + delta;
				P(i,j) = P(i,j) + 2*div;
				P(i+1,j) = P(i+1,j) - div/4;
				P(i,j+1) = P(i,j+1) - div/4;
			elseif i== 1 && j==N
				U(i+1,j) = U(i+1,j) + delta;
				V(j,i) = V(j,i) - delta;
				P(i,j) = P(i,j) + 2*div;
				P(i+1,j) = P(i+1,j) - div/4;
				P(i,j-1) = P(i,j-1) - div/4;
			elseif i ==N && j==1
				U(i,j) = U(i,j) - delta;
				V(j+1,i) = V(j+1,i) + delta;
				P(i,j) = P(i,j) + 2*div;
				P(i-1,j) = P(i-1,j) - div/4;
				P(i,j+1) = P(i,j+1) - div/4;
			else
				U(i,j) = U(i,j) - delta;
				V(j,i) = V(j,i) - delta;
				P(i,j) = P(i,j) + 2*div;
				P(i-1,j) = P(i-1,j) - div/4;
				P(i,j-1) = P(i,j-1) - div/4;
			end
		end
	end
	U = U(2:N,:);
	V = V(2:N,:);


