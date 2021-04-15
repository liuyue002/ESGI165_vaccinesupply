function   DX = odefun_SRI(t, X, xi, M, N, beta, gamma)
unvaccinated=N-X(3*M+1:4*M);
unvaccinated(unvaccinated==0)=1; % avoid division by 0
DX = [-beta*X(M+1:2*M).*X(1:M)./N - xi(t).*X(1:M)./unvaccinated ; 
       beta*X(M+1:2*M).*X(1:M)./N - gamma*X(M+1:2*M) - xi(t).*X(M+1:2*M)./unvaccinated ;
       gamma*X(M+1:2*M) - xi(t).*X(2*M+1:3*M)./unvaccinated ;
       xi(t)];