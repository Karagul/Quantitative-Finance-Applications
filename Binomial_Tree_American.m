
%% Binomial Trees American Style

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
American = true; 
% InnerValue =@(S,K) max(S-K,0);
InnerValue =@(S,K) max(K-S,0); % put option 

%% prepare binomail tree
dt = T/M;                              % delta t, time increments/ length of 1 subperiod in years 
d = exp(-sigma*sqrt(dt));
u = 1/d;
p = (exp(r*dt)-d)/ (u-d);   % p = 0.4941 likelihood of upmove$


%offset for indexing with j and i

f7 = 1; 

%% stock price at time T

j = 0:M; 
ST = S0 * (u.^j .* d.^(M-j));
% STalt = S0 * u.^(-M:2:M);                   % for case u = 1/d, other
                                             %implementation

plot(ST, '.') % all of these are possible

%binomial expansion 

%  for j=0:M
%     bn(j) = nchoosek(M,j);   
% end


for j=0:M
        bn(j+f7) = nchoosek(M,j); 
end

j = 0:M; % redefine j as a vector
prob = bn .* p.^j .* (1-p).^(M-j);

plot(ST,prob)

% sum(prob)-1                   % should ad up to one

% ST * prob' .                  
 

%% option pricing (American Style)


% Stock prices 

S_t = nan(M+f7, M+f7);

S_t(0+f7,0+f7) = S0;            %root of tree at (0,0)
for i = 1:M
     j = 0:(i-1);
     S_t(i +f7, j +f7) = S_t(i-1 +f7, j   +f7) *d;
     S_t(i +f7, i +f7) = S_t(i-1 +f7, i-1 +f7) *u;
    
end
S_t

plot(0:M, S_t, '*k'); % increase M 

semilogy(0:M, S_t, '*k'); % creates semilog plot log on y  

%% pricing the option

opt_t = nan(size(S_t)); 

% at maturity 
opt_t(M +f7, :) = InnerValue(S_t(M +f7, :), K); 

% towards t=0
for i = (M-1):-1:0
    j = 0:i; % safes another loop 
    opt_t( i +f7, j +f7) = (  p    * opt_t(i+1 +f7, j+1 +f7) ...      % up   % expacted value
                           + (1-p) * opt_t(i+1 +f7, j   +f7)  ) ...   % down
                           * exp(-r*dt);           
               if American 
                   vie = InnerValue(S_t(i +f7, j +f7), K); % value if we exercise right now
                   opt_t(i +f7, j +f7) = max(opt_t(i +f7, j +f7), vie); 
               end
               
         
end

% results: european / "grow the tree (Eur/Am)"
disp([option_0, opt_t(0 +f7, 0 +f7)]) % if the same then pre-mature exercise has no value, for put option has extra value, change to put then should be different

%should be same
opt_t(1,1);
option_0;



