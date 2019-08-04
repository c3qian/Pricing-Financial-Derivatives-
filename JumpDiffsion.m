
% This code is used to simulate jump diffusion
% assume no more than 1 jump happens in [t, t+dt]

randn(’state’,100);
rand(’state’, 10);

T = 1.00;% expiry time
sigma = 0.25; % volitility
mu = .10; %drift
S_init = 100;


% assume jump size: log normal distribution

lambda = .1 ;%jump size arrival rate
jump_vol = .40; %standard deviation of jump size
jump_mean = -.9; % mean of the jump size

N_sim = 100000; % total num of simulations
N = 250; % num of steps
delt = T/N; % delta t 

% compensated  drift : E[j-1]
kappa = exp(.5*jump_vol*jump_vol + jump_mean) - 1.;

% compensated drift for X = log(S)

drift = (mu - sigma*sigma/2.0 - lambda*kappa);

% X = log(S)

X_old(1:N_sim,1) = log(S_init); X_new(1:N_sim,1) = zeros(N_sim, 1);

jump_chek = zeros(N_sim, 1); jump_size = zeros(N_sim, 1); jump_mask = zeros( N_sim, 1);

for i=1:N % timestep loop

	jump_chek(:,1) = rand(N_sim,1); jump_mask(:,1) = ( jump_chek(:,1) <= lambda*delt);

	jump_size(:,1) = jump_mean + jump_vol*randn(N_sim,1); jump_size = jump_size.*jump_mask;

	X_new(:,1) = X_old(:,1) + drift*delt +sigma*sqrt(delt)*randn(N_sim,1) +jump_size(:,1);

	X_old(:,1) = X_new(:,1);

end % timestep loop

% visualization
S(:,1) = exp( X_new(:,1) );

n_bin = 200; 
hist(S, n_bin);

stndrd_dev = std(S); 
disp(sprintf(’standard deviation:%.5g\n’,stndrd_dev));

mean_S = mean(S); 
disp(sprintf(’mean:%.5g\n’,mean_S))
