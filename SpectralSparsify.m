function H = SpectralSparsify(G,epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs: -G           similarity matrix of the original graph of size [nxn]
%%%         -epsilon     sparsification parameter
%%% Output: -H           similarity matrix of the sparsified graph [nxn]
%%%%%%%%%%%%%%%%%%%%%%%%

n = size(G,1);
N = ceil(8*n*log(n)/epsilon^2); %%% number of edge samples
nonzero = find(triu(G));
[s,d] = find(triu(G));
elist = [s d];
w = G(nonzero);
R = EffectiveResistances(elist,elist,w,1e-5,1,'slm'); %%% you can dowslonad this function and its dependence from http://www.cs.cmu.edu/~jkoutis/SpectralAlgorithms.htm
p = w.* R; p = p/sum(p); %%% probability of the edges sampling


H = zeros(n,n);
f = randsample(size(p,1),N,true,p); %%% sampling N edges based on p
unique_f =  unique(f);
count_f = hist(f,unique_f);
k = nonzero(unique_f);

notrep = find(count_f == 1);
H(k(notrep)) = w(unique_f(notrep))./(N*p(unique_f(notrep)));

rep = find(count_f>1);
for i = 1:size(rep,2)
    H(k(rep(i))) = count_f(rep(i))*w(unique_f(rep(i)))/(N*p(unique_f(rep(i))));
end

H = H + H';

end

