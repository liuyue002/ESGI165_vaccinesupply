M=3; % number of countries
gamma=1/14;
%beta=(eye(M)*1.6 + 0.03)*gamma;
beta=[1.6,0,0;
      0.03,1.6,0;
      0.03,0,1.6]*gamma;
% beta=(eye(M)*1.6 + 0.03*2*rand(M))*gamma;

%%
%xi=[0.00;0.0;0.0];
N=[6;3.6;3];
%xi=[0.0;0.0028;0.0033]; %constaint: sum(xi)*N=0.02
%xi=@(t) [0.01;0.01;0.01];

Sinit=[0.9; 1; 0.5].*N;
Iinit=[0.1; 0; 0].*N;
Rinit=[0;0;0].*N;
Vinit=[0;0;0.5].*N;

tmax = 200;


K = sum(N)/200; % Total number of vaccines

k=0.05*N;
m1=1;
m2=2;
C1 = @(xi) cost_SIR(xi, M, beta, gamma, N, Sinit, Iinit, Rinit, Vinit, tmax,m1,m2,k);

Aeq = ones(1,M);

xi0 = (ones(M,1)+10^(-2)*randn(M,1));
xi0 = K/(sum(xi0)) * xi0;

A=[]; b=[]; 
lb=zeros(M,1);


options = optimoptions('fmincon');


options.Display = 'off';
options.UseParallel = false;

[x,fval, exitflag] = fmincon(C1,xi0,A,b,Aeq,K, lb, b, b, options);