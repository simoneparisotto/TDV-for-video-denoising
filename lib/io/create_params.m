function params = create_params(alpha,sigma,fps,ST_SIGMAS,ST_RHOS,ETA,maxiters,verbose,tolerance)

stdev    = sigma/255;
variance = stdev^2;

%mu  = 10;  % weight of the regulariser
%eta = (mu*variance)^-1;
eta = ETA;

% common parameters
params.varying      = 1;
params.lambda       = alpha;
params.order        = [1 2 3];
params.stdev        = stdev;
params.variance     = variance;
params.eta          = eta;
params.acceleration = 1;
params.verbose_text = verbose;

params.sigma        = ST_SIGMAS;
params.rho          = ST_RHOS;
params.time_gain    = 1;
params.FrameRate    = fps;
params.tolerance    = tolerance;
params.niter        = maxiters;