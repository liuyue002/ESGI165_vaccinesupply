function [C,I] = cost_SIR(xi, M, beta, gamma, N, Sinit, Iinit, Rinit, Vinit, tmax)
   
           
    % Generate function xi -- change as needed!
    
    xi = @(t) xi;
    
    tspan=[0,tmax];
    init=[Sinit;Iinit;Rinit;Vinit];

    options = odeset('RelTol',1e-12,'AbsTol',1e-13,'Events',@(t,X) event_negative(t,X,M));

    IsItEnd = 0;
    t = [];
    X = [];
    while IsItEnd == 0
        % Run the ODE solver
        [tt,XX,te,ye,ie]=ode45(@(tt,XX) odefun_SRI(tt, XX, xi, M, N, beta, gamma),tspan,init,options);

        % Determine why the solver halted
        if isempty(ie) % Run out of time
            
            t = [t; tt];
            X = [X; XX];
            IsItEnd = 1;
        else % One of Susceptible populations reached 0
            for m = ie % Cancel all transmitions and vaccinations in that population
                
                t = [t; tt];
                X = [X; XX];

                beta(:,m) = zeros(M,1);
                beta(m,:) = zeros(1,M);

                % Redefine function xi
                mmultip = ones(M,1);
                mmultip(m) = 0;
                xi = @(t) mmultip.*xi(t);

                % Redefine ODE problem
                tspan = [t(end), tmax];
                Sinit = X(end,1:M);
                Iinit = X(end,M+1:2*M);
                Rinit = X(end,2*M+1:3*M);
                Vinit = X(end,3*M+1:4*M);

                clear tt XX
                % Ensure no more susceptibles or recovered can be found
                Sinit(m) = 0;
                Rinit(m) = 0;
                Iinit(m) = 0;

                init=[Sinit,Iinit,Rinit,Vinit];

            end
        end
    end
    
    I = trapz(t,X(:,M+1:2*M));
    C = sum(I);
end
