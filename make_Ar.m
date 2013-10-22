function A=make_Ar(r)

% This program accompanies the paper Rank fooling set size by 
% Aya Hamed and Troy Lee

% Given input r, the program outputs a matrix A of rank r and size
% r+1 choose 2 that is a Hadamard factorization of the identity, i.e.
% A.*A'==eye(nchoosek(r+1,2));

R=[0:2*r-1];
C=[-1:r-1];
P=binom(R,C);
v=[0:r-1];
v=(-1).^v;
P=[[1,v];P];

% P is an extended Pascal's triangle, with nonstandard def -1 choose -1=1
% we use P to make the F_k matrices with toeplitz command

A=[];
for i=0:r-1
	M=[];
	for j=0:r-1
		k=j-i;
		if k >0
			M=[M;toeplitz(zeros(1,r-j),P(k+1,1:r-i))];
		else
			k=-k;
			row=v(1:r-i).*P(k+1:k+r-i,k+1)';
			col=(-1)^k*v(1:r-j).*P(1:r-j,k+2)';
			M=[M;toeplitz(col,row)];
		end;
	end;
	A=[A,M];
end
