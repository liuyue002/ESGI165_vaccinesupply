function [C,costpercountry] = cost_SIR(xi_params, M, beta, gamma, N, Sinit, Iinit, Rinit, Vinit, tmax,m1,m2,k,costpervaccine,make_plot)


% Generate function xi -- change as needed!

% assume xi is piecewise constant, with 3 pieces, cutoff at t1, t2
% column 1 of xi_params are the xi's for the first time range, and so on
%t1=50;
%t2=100;
%xi_params = reshape(xi_params,M,3);
%xi = @(t) xi_params(:,1) + heaviside(t-t1)*(xi_params(:,2)-xi_params(:,1)) + heaviside(t-t2)*(xi_params(:,3)-xi_params(:,2));

%xi_params(3) is the amount that country 2 is given to country 1 for free
%country 2 always produce maximum amount 0.1
xi=@(t) [xi_params(1)+xi_params(3);xi_params(2)-xi_params(3)];
xi_original = xi;

tspan=[0,tmax];
init=[Sinit;Iinit;Rinit;Vinit];

options = odeset('RelTol',1e-12,'AbsTol',1e-13,'Events',@(t,X) event_negative(t,X,M));

IsItEnd = 0;
t = [];
X = [];
tHalt = zeros(M,1);
tHaltIndex = zeros(M,1);
while IsItEnd == 0
    % Run the ODE solver
    [tt,XX,te,ye,ie]=ode45(@(tt,XX) odefun_SRI(tt, XX, xi, M, N, beta, gamma),tspan,init,options);
    
    % Determine why the solver halted
    if isempty(ie) % Run out of time
        
        t = [t; tt];
        X = [X; XX];
        IsItEnd = 1;
    else % One of Susceptible populations reached 0
        t = [t; tt];
        X = [X; XX];
        tspan = [t(end), tmax];
        if t(end)>=tmax
            IsItEnd = 1;
        end
        Sinit = X(end,1:M);
        Iinit = X(end,M+1:2*M);
        Rinit = X(end,2*M+1:3*M);
        Vinit = X(end,3*M+1:4*M);
        for m = ie' % Cancel all transmitions and vaccinations in that population
            tHalt(m) = tt(end);
            tHaltIndex(m) = size(t,1);
            beta(:,m) = zeros(M,1);
            beta(m,:) = zeros(1,M);
            
            % Redefine function xi
            mmultip = ones(M,1);
            mmultip(m) = 0;
            xi = @(t) mmultip.*xi(t);
            
            
            % Ensure no more susceptibles or recovered can be found
            Sinit(m) = 0;
            Rinit(m) = 0;
            Iinit(m) = 0;
            
           
            
        end
        clear tt XX
        init=[Sinit,Iinit,Rinit,Vinit];
    end
end

% death rate, k=hospital capacity, m=additional penalty past k
drate=@(I,k)(m1*I + (m2-m1)*(I>k) );
Itotal=zeros(M,1);
for i=1:M
    Itotal(i)=trapz(t,drate(X(:,M+i),k(i)));
end
%vaccinecost=costpervaccine.*(X(end,3*M+1:4*M)');
vaccinecost=costpervaccine.*[xi_params(1);xi_params(2)];
%productioncost=zeros(M,1);

costpercountry=Itotal+vaccinecost;

C = sum(costpercountry);

if make_plot
    fig=figure('Position',[121 346 1734 439]);
    for i=1:M
        subplot(1,M,i);
        hold on
        plot(t,X(:,i)./N(i));
        plot(t,X(:,M+i)./N(i));
        plot(t,X(:,2*M+i)./N(i));
        plot(t,X(:,3*M+i)./N(i));
%         boundaryline=xline(t1);
%         boundaryline.Alpha=0.3;
%         boundaryline=xline(t2);
%         boundaryline.Alpha=0.3;
        horizline=yline(k(i)/N(i));
        horizline.Alpha=0.3;
        xlabel('t');
        axis([0 tmax 0 1]);
        legend('S','I','R','V');
        %fprintf('integral of I1: %.3f\n', Itotal(1));
        title(['country ',num2str(i),', Itotal=',num2str(Itotal(i),'%.3f'),' ,vaccinecost=',num2str(vaccinecost(i))]);
        hold off
    end
end

end
