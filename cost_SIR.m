function [C,Itotal] = cost_SIR(xi_params, M, beta, gamma, N, Sinit, Iinit, Rinit, Vinit, tmax,m1,m2,k,make_plot)


% Generate function xi -- change as needed!

% assume xi is piecewise constant, with 3 pieces, cutoff at t1, t2
% column 1 of xi_params are the xi's for the first time range, and so on
t1=30;
t2=60;
xi_params = reshape(xi_params,M,3);
xi = @(t) xi_params(:,1) + heaviside(t-t1)*(xi_params(:,2)-xi_params(:,1)) + heaviside(t-t2)*(xi_params(:,3)-xi_params(:,2));

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

% death rate, k=hospital capacity, m=additional penalty past k
drate=@(I,k)(m1*I + (m2-m1)*(I>k) );
Itotal=zeros(M,1);
for i=1:M
    Itotal(i)=trapz(t,drate(X(:,M+i),k(i)));
end
C = sum(Itotal);

if make_plot
    fig=figure('Position',[121 346 1734 439]);
    sfig1=subplot(1,3,1);
    hold on
    plot(t,X(:,1)./N(1));
    plot(t,X(:,4)./N(1));
    plot(t,X(:,7)./N(1));
    plot(t,X(:,10)./N(1));
    boundaryline=xline(t1);
    boundaryline.Alpha=0.3;
    boundaryline=xline(t2);
    boundaryline.Alpha=0.3;
    xlabel('t');
    axis([0 tmax 0 1]);
    legend('S1','I1','R1','V1');
    intI=trapz(t,X(:,4));
    fprintf('integral of I1: %.3f\n', intI);
    title(['country 1, intI=',num2str(intI,'%.3f')]);
    
    sfig2=subplot(1,3,2);
    hold on
    plot(t,X(:,2)./N(2));
    plot(t,X(:,5)./N(2));
    plot(t,X(:,8)./N(2));
    plot(t,X(:,11)./N(2));
    boundaryline=xline(t1);
    boundaryline.Alpha=0.3;
    boundaryline=xline(t2);
    boundaryline.Alpha=0.3;
    xlabel('t');
    axis([0 tmax 0 1]);
    legend('S2','I2','R2','V2');
    intI=trapz(t,X(:,5));
    fprintf('integral of I2: %.3f\n', intI);
    title(['country 2, intI=',num2str(intI,'%.3f')]);
    
    sfig3=subplot(1,3,3);
    hold on
    plot(t,X(:,3)./N(3));
    plot(t,X(:,6)./N(3));
    plot(t,X(:,9)./N(3));
    plot(t,X(:,12)./N(3));
    boundaryline=xline(t1);
    boundaryline.Alpha=0.3;
    boundaryline=xline(t2);
    boundaryline.Alpha=0.3;
    xlabel('t');
    title('country 3');
    axis([0 tmax 0 1]);
    legend('S3','I3','R3','V3');
    intI=trapz(t,X(:,6));
    fprintf('integral of I3: %.3f\n', intI);
    title(['country 3, intI=',num2str(intI,'%.3f')]);
end

end
