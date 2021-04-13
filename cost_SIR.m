function [C,I1,I2,I3] = cost_SIR(xi, M, beta, gamma, N, Sinit, Iinit, Rinit, Vinit, tmax)
    odefun=@(t,X)[-beta*X(M+1:2*M).*X(1:M)./N - xi.*X(1:M)./(N-X(3*M+1:4*M)) ; 
               beta*X(M+1:2*M).*X(1:M)./N - gamma*X(M+1:2*M) - xi.*X(M+1:2*M)./(N-X(3*M+1:4*M)) ;
               gamma*X(M+1:2*M) - xi.*X(2*M+1:3*M)./(N-X(3*M+1:4*M)) ;
               xi];
           
    tspan=[0,tmax];
    init=[Sinit;Iinit;Rinit;Vinit];

    options = odeset('RelTol',1e-10,'AbsTol',1e-7);

    [t,X]=ode45(odefun,tspan,init,options);
    I1 = trapz(t,X(:,4));
    I2 = trapz(t,X(:,5));
    I3 = trapz(t,X(:,6));
    C =  I1+I2+I3;
end
