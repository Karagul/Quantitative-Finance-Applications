%% Binomial Trees 
clear 

%% settings

%Stock
S0 = 100; 
r = 0.02;                % return p. a.
sigma = 0.30;            % vola p. a.

% time horizon
T = 1;                   % time to maturity (in years)
M = 50;                  % number of subperiods from 0 to T (in trading days)

% option
K = 100; 
InnerValue =@(S,K) max(S-K,0);


%% prepare binomail tree
dt = T/M;                              % delta t, time increments/ length of 1 subperiod in years 
d = exp(-sigma*sqrt(dt));
u = 1/d;
p = (exp(r*dt)-d)/ (u-d);   % p = 0.4941 likelihood of up move


%offset for indexing with j and i

f7 = 1; 

%% stock price at time T

j = 0:M; 
ST = S0 * (u.^j .* d.^(M-j));
% STalt = S0 * u.^(-M:2:M);                   % for case u = 1/d, other
                                             %implementation

% plot(ST, '.') % all of these are possible

% binomial expansion 

%  for j=0:M
%     bn(j) = nchoosek(M,j);    
% end


for j=0:M
        bn(j+f7) = nchoosek(M,j); 
end

j = 0:M; % redefine j as a vector
prob = bn .* p.^j .* (1-p).^(M-j); 

plot(ST,prob)

% sum(prob)-1            

% ST * prob' .                 
% S0 * exp(r*T)
 

%% option pricing

option_T = InnerValue(ST,K);

%plot(ST, option, '.');

option_0 = (option_T * prob' ) * exp(-r*T);
option_0                                         













