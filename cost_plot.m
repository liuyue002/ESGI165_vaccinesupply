M=2; % number of countries
gamma=1/14;
%beta=(eye(M)*1.6 + 0.03)*gamma;
beta=[1.6,0.3;
      0.3,1.6]*gamma;
%  beta=(eye(M)*1.6 + 0.03*2*rand(M))*gamma;

%%
%xi=[0.00;0.0;0.0];
N=[6;3];
%xi=[0.0;0.0028;0.0033]; %constaint: sum(xi)*N=0.02
%xi=@(t) [0.01;0.01;0.01];

Sinit=[0.9; 1].*N;
Iinit=[0.1; 0].*N;
Rinit=[0;0].*N;
Vinit=[0;0].*N;

tmax = 300;


K = 0.08; % Total number of vaccines

numpts=30;
xi1=linspace(0,K,numpts);
xi2=linspace(0,K,numpts);
[Xi1,Xi2]=meshgrid(xi1,xi2);

k=0.05*N;
m1=1;
m2=2;
%costpervaccine=[4;2];
costpervaccine=[8;4];
C1 = @(xi1,xi2) cost_SIR([xi1;xi2], M, beta, gamma, N, Sinit, Iinit, Rinit, Vinit, tmax,m1,m2,k,costpervaccine,0);
[cost,costpercountry]=arrayfun(C1,Xi1,Xi2,'UniformOutput',false);
cost=cell2mat(cost);
cost1=zeros(size(costpercountry));
cost2=zeros(size(costpercountry));
for i=1:numpts
    for j=1:numpts
        cost1(i,j)=costpercountry{i,j}(1);
        cost2(i,j)=costpercountry{i,j}(2);
    end
end
%%
fig=figure;

subplot(1,3,1);
surf(Xi1,Xi2,cost);
xlabel('\xi_1');
ylabel('\xi_2');
zlabel('cost');
title('total cost');

subplot(1,3,2);
surf(Xi1,Xi2,cost1);
xlabel('\xi_1');
ylabel('\xi_2');
zlabel('cost');
title('cost of country 1');

subplot(1,3,3);
surf(Xi1,Xi2,cost2);
xlabel('\xi_1');
ylabel('\xi_2');
zlabel('cost');
title('cost of country 2');