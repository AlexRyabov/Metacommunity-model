function Y = MDS(D, nl);
N=length(D);
subB = -.5*(D.^2 - sum(D'.^2)'*ones(1,N)/nl - ones(N,1)*sum(D.^2)/N+sum(sum(D.^2))/(N*nl));
[alpha,beta] = eigs(subB'*subB, nl, 'LR');
val = beta.^(1/2);
vec = subB*alpha*inv(val);
h = real(diag(val));
[foo,sorth] = sort(h);  
sorth = sorth(end:-1:1);
val = real(diag(val(sorth,sorth)));
vec = vec(:,sorth);
Y = real(vec(:,1:nl).*(ones(N,1)*sqrt(val(1:nl))'))';





