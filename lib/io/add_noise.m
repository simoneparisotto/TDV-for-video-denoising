function u = add_noise(uorig,sigma)

%% ADD NOISE in [0,1]
% sigma is the standard deviation in [0,255]
% stdev is the standard deviation in [0,1]
% var   is the variance in [0,1]
rng('default')
seed   = 0; rng(seed);
R      = randn(size(uorig));
% Add noise from normal distribution with mean 0 and standard deviation sigma/255.
u      = uorig + (sigma/255)*R;

end