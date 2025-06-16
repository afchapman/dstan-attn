% this codes generates model output used for contrast-dependent suppression
% analyses and was run on the BU SCC

% data output is available on our OSF repository: https://osf.io/qy9pa/

ncores = str2num(getenv("NSLOTS"));
pool = parpool(ncores);

addpath('model');

%% model setup
opt = [];
modelClass = [];
rsoa = 2001; % SOA = 250 ms (see runModel)
rseq = 2; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.aAI = 0;
opt.aAV = 0;
opt.sigma1 = 0.02;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

opt.tauE1 = 400;
opt.tauS1 = 100;

opt.dt = 2;
opt.T = 2.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.stimDur = 100;

opt.x = [-1 1];
opt.nx = 2;

opt.stimMode = 'surr_supp';

%% generate model output for each parameter combination
contrC = [0 logspace(log10(0.05),log10(1),9)]; % center
contrS = [0 .12 .25 .5 1]; % surround

paramList = combvec(contrC,contrS);

r1_surr = nan(length(paramList),24,opt.nt);
parfor ii=1:length(paramList)
    opt2 = opt;
    opt2.stimContrasts = paramList(:,ii);
    [~,p,~] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    r1_surr(ii,:,:) = p.r1;
end

save('output/surroundSupp.mat','r1_surr');
