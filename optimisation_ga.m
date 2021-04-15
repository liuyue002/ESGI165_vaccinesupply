M=3; % number of countries
gamma=1/14;
%beta=(eye(M)*1.6 + 0.03)*gamma;
beta=[1.6,0,0;
      0.03,1.6,0;
      0.03,0,1.6]*gamma;
%  beta=(eye(M)*1.6 + 0.03*2*rand(M))*gamma;

%%
%xi=[0.00;0.0;0.0];
N=[6;3.6;3];
%xi=[0.0;0.0028;0.0033]; %constaint: sum(xi)*N=0.02
%xi=@(t) [0.01;0.01;0.01];

Sinit=[0.9; 1; 1].*N;
Iinit=[0.1; 0; 0].*N;
Rinit=[0;0;0].*N;
Vinit=[0;0;0].*N;

tmax = 200;


K = sum(N)/200; % Total number of vaccines

k=0.05*N;
m1=1;
m2=2;
costpervaccine=[4;2;2];
C1 = @(xi_params) cost_SIR(xi_params, M, beta, gamma, N, Sinit, Iinit, Rinit, Vinit, tmax,m1,m2,k,costpervaccine,0);

A = zeros(3,M*3);
for i=1:3
    A(i,(i-1)*M+1:i*M)=1;
end

% assume there are 3 time ranges, cutoff at t1, t2
num_timerange=3;
xi0 = (ones(M,num_timerange)+10^(-2)*randn(M,num_timerange));
xi0 = (K./sum(xi0,1)).*xi0;
xi0 = reshape(xi0,M*3,1);

Aeq=[];
lb=zeros(M*num_timerange,1);


options = optimoptions('ga');


options.Display = 'off';
options.UseParallel = true;
options.InitialPopulationMatrix = xi0';

[x,fval, exitflag] = ga(C1,length(xi0),Aeq,ones(num_timerange,1)*K, A,b, lb, ones(size(xi0))*K, b, options);

%%
optimal_xi=x;
cost_SIR(optimal_xi, M, beta, gamma, N, Sinit, Iinit, Rinit, Vinit, tmax,m1,m2,k,1);