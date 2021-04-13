function   DX = odefun_SRI(t, X, xi, M, N, beta, gamma)

DX = [-beta*X(M+1:2*M).*X(1:M)./N - xi(t).*X(1:M)./(N-X(3*M+1:4*M)) ; 
       beta*X(M+1:2*M).*X(1:M)./N - gamma*X(M+1:2*M) - xi(t).*X(M+1:2*M)./(N-X(3*M+1:4*M)) ;
       gamma*X(M+1:2*M) - xi(t).*X(2*M+1:3*M)./(N-X(3*M+1:4*M)) ;
       xi(t)];