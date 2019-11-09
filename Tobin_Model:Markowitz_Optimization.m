%% Markowitz
clear; clf; clc; 
% loading the empirical data



rF = 0.00; 


S = dlmread('DJI_Dow_Jones_Industrial_Average.csv', ',', 1, 0);

disp(datestr(S([1, end],1))); 

S(:,1:2) = [];

a = [1:3];  % assets 
N = numel(a); 
%% estimate parameters 

rEmp =  diff(log(S(:, a)));  % compute log returns 

% parameters (per day)
E = mean(rEmp); 
SD = std(rEmp); 
V = cov(rEmp); 

plot(SD, E, '.', 'markersize', 30)



%% random pf 


nExp = 10000; 
x = rand(N, nExp); % need to sum up to two sum(x/sum(x))                
x = x ./ sum(x);
rP = E *x;  % expected return of portfolie
sP = sqrt(diag(x' * V * x)); % volatility 

% SD = 0.0154    0.0130    0.0151 %sqrt(x' * V * x) =0.0112 % pf has less
% risk then indiv assets

hold on
plot(sP,rP, '*');
hold off
xlim([.01, .017])


% cardinality constraint 

nExp = 1000; 
EPew = nan(nExp,N);
VPew = nan(nExp,N); 

for n = 1:N % Number of different assets in the pf (cardinality)
    x = ones(n,1)/n;

    for e = 1:nExp
         i = randperm(N,n); 
         Epew(e,n) = E(i)*x; % compute expected return...
         VPew(e,n) = sqrt(x' * V(i,i) * x);
    
    end
end


% subplot(2,1,1)
% boxplot(EPew); xlabel('cardinality'); ylabel('expected return')
% subplot(2,1,2)
% boxplot(VPew); xlabel('cardinality'); ylabel('volatility') % beneficial to have large cardinality, even if you pick the worst combination you are better of then single 
% 


%% Markowitz optimization 



H = 2*V;
f = zeros(N,1);
A = [];
b = [];
Aeq = ones(1,N); % sum up all weights
beq = 1; %update this one
ub = ones(N,1)*1000; % no upper bound, assign very large value
lb = zeros(N,1); % lower bound


xMVP = quadprog(H, f, A, b, Aeq, beq, lb, ub) % find the optimal weights

rMVP = E*xMVP; % compute MVP's expected return 
sMVP = sqrt(xMVP' * V * xMVP);
hold on 
plot(sMVP, rMVP, '*k', 'markersize', 20); 
hold off


% efficient line 

nSteps = 51; 
rTarget = linspace(rMVP, max(E), nSteps);
xM= nan(N, nSteps);

for i =1:numel(rTarget) % for each target return candidate 
    Aeq = [ones(1,N);  % 1st eq-contraint: sum up all the weights
        E ];         % 2end eq-constraint: for target return 
    beq = [ 1;
        rTarget(i)];
   
    xM(:,i) = quadprog(H, f, A, b, Aeq, beq, lb, ub);
end
rPM = E*xM; 
sPM = sqrt(diag(xM' * V * xM ) ); 

hold on 
plot(sPM, rPM, 'k', 'linewidth', 3); 
hold off



%% Tobin Model 
 


% fmincon =finction minmization constrain


% objective: 
negTheta =@(x) -(  (E*x - rF) / sqrt(x'*V*x)   );
% weights have to sum up to 1
Aeq = ones(1,N); 
beq =1; 
% lower bound for x_i: (here: 0)
A = -eye(N); 
b = ones(N,1)* 0;
% initial solution: equal weights
x0 = ones(N,1) / N; 

% find tangency portfolio 
xT = fmincon(negTheta, x0, A, b, Aeq, beq)

rT = E * xT; 
sT = sqrt(xT' * V * xT); 

hold on 
    plot(sT, rT, '.k', 'markersize', 30);
    plot([0, sT], [rF, rT], 'k');
 hold off




 
 
 
 