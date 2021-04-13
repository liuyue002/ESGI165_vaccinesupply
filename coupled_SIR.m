
M=3; % number of countries
gamma=1/14;
%beta=(eye(M)*1.6 + 0.03)*gamma;
% beta=[1.6,0,0;
%       0.03,1.6,0;
%       0.03,0,1.6]*gamma;
%beta=(eye(M)*1.6 + 0.03*2*rand(M))*gamma;

%%
%xi=[0.00;0.0;0.0];
N=[6;3.6;3];
%xi=[0.0;0.0028;0.0033]; %constaint: sum(xi)*N=0.02
xi=@(t) [0.01;0.01;0.01];
% S: X(1:M), I: X(M+1:2*M), R: X(2*M+1:3*M), V: X(3*M+1:4*M)
odefun=@(t,X)[-beta*X(M+1:2*M).*X(1:M)./N - xi(t).*X(1:M)./(N-X(3*M+1:4*M)) ; 
               beta*X(M+1:2*M).*X(1:M)./N - gamma*X(M+1:2*M) - xi(t).*X(M+1:2*M)./(N-X(3*M+1:4*M)) ;
               gamma*X(M+1:2*M) - xi(t).*X(2*M+1:3*M)./(N-X(3*M+1:4*M)) ;
               xi(t)];
           
tspan=[0,500];
Sinit=[0.9; 1; 1].*N;
Iinit=[0.1; 0; 0].*N;
Rinit=[0;0;0].*N;
Vinit=[0;0;0].*N;
init=[Sinit;Iinit;Rinit;Vinit];

options = odeset('RelTol',1e-12,'AbsTol',1e-13,'Events',@(t,X) event_negative(t,X,M));

[t,X]=ode45(odefun,tspan,init,options);

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
legend('S3','I3','R3','V3');
fprintf('integral of I3: %.3f\n', trapz(t,X(:,6)));