function AI = AmariIndex(G)

G = abs(G);
N = size(G,1);

AI = 1/(2*N*(N-1))*sum(sum(abs(G)./repmat(max(abs(G),[],1),N,1),1)-1) ...
    + 1/(2*N*(N-1))*sum(sum(abs(G)./repmat(max(abs(G),[],2),1,N),2)-1);

