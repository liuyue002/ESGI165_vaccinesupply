
M=3; % number of countries
gamma=1/14;
%beta=(eye(M)*1.6 + 0.03)*gamma;
% beta=[1.6,0,0;
%       0.03,1.6,0;
%       0.03,0,1.6]*gamma;
beta=(eye(M)*1.6 + 0.03*2*rand(M))*gamma;
beta1 = beta;
%%
%xi=[0.00;0.0;0.0];
N=[6;3.6;3];
%xi=[0.0;0.0028;0.0033]; %constaint: sum(xi)*N=0.02
xi=@(t) [0.02;0.01;0.01];
% S: X(1:M), I: X(M+1:2*M), R: X(2*M+1:3*M), V: X(3*M+1:4*M)


tmax = 700;
tspan=[0,tmax];
Sinit=[0.9; 1; 1].*N;
Iinit=[0.1; 0; 0].*N;
Rinit=[0;0;0].*N;
Vinit=[0;0;0].*N;
init=[Sinit;Iinit;Rinit;Vinit];

options = odeset('RelTol',1e-12,'AbsTol',1e-13,'Events',@(t,X) event_negative(t,X,M));

IsItEnd = 0;
t = [];
X = [];
z = 0;
while IsItEnd == 0
    % Run the ODE solver
    [tt,XX,te,ye,ie]=ode45(@(tt,XX) odefun_SRI(tt, XX, xi, M, N, beta, gamma),tspan,init,options);
    ie = ie';
    % Determine why the solver halted
    if isempty(ie) % Run out of time
        
        t = [t; tt];
        X = [X; XX];
        IsItEnd = 1;
    else % One of Susceptible populations reached 0
        t = [t; tt];
        X = [X; XX];
        tspan = [t(end), tmax];
            Sinit = X(end,1:M);
            Iinit = X(end,M+1:2*M);
            Rinit = X(end,2*M+1:3*M);
            Vinit = X(end,3*M+1:4*M);
        for m = ie % Cancel all transmitions and vaccinations in that population
            m
            
            
            beta(:,m) = zeros(M,1);
            beta(m,:) = zeros(1,M);
            
            % Redefine function xi
            mmultip = ones(M,1);
            mmultip(m) = 0;
            xi = @(t) mmultip.*xi(t);
            
            % Redefine ODE problem
            
            
            
            % Ensure no more susceptibles or recovered can be found
            Sinit(m) = 0;
            Rinit(m) = 0;
            Iinit(m) = 0;
            z = z+1;
            
            
        end
        init=[Sinit,Iinit,Rinit,Vinit];
        clear tt XX
    end
    if z == M
        IsItEnd = 1;
    end
end
%%
fig=figure('Position',[121 346 1734 439]);
sfig1=subplot(1,3,1);
hold on
plot(t,X(:,1)./N(1));
plot(t,X(:,4)./N(1));
plot(t,X(:,7)./N(1));
plot(t,X(:,10)./N(1));
xlabel('t');
title('country 1');
axis([0 tmax 0 1]);
legend('S1','I1','R1','V1');
fprintf('integral of I1: %.3f\n', trapz(t,X(:,4)));

sfig2=subplot(1,3,2);
hold on
plot(t,X(:,2)./N(2));
plot(t,X(:,5)./N(2));
plot(t,X(:,8)./N(2));
plot(t,X(:,11)./N(2));
xlabel('t');
title('country 2');
axis([0 tmax 0 1]);
legend('S2','I2','R2','V2');
fprintf('integral of I2: %.3f\n', trapz(t,X(:,5)));

sfig3=subplot(1,3,3);
hold on
plot(t,X(:,3)./N(3));
plot(t,X(:,6)./N(3));
plot(t,X(:,9)./N(3));
plot(t,X(:,12)./N(3));
xlabel('t');
title('country 3');
axis([0 tmax 0 1]);
legend('S3','I3','R3','V3');
fprintf('integral of I3: %.3f\n', trapz(t,X(:,6)));