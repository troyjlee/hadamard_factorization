function M=binom(n,k)

% This function is a matricized version of nchoosek
% Given row vectors n,k the output is a matrix M whose (i,j) is 
% n(i) choose k(j)

s=size(n,2);
t=size(k,2);

N=n'*ones(1,t);
K=ones(s,1)*k;
M=N-K+1;
M=gamma(M);
M=M.*gamma(K+1);
M=gamma(N+1)./M;
