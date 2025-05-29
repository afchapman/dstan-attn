% this codes generates model output used for response adaptation analyses
% with orthogonal stimuli in a sequence and was run on the BU SCC

% data output is available on our OSF repository: https://osf.io/qy9pa/

ncores = str2num(getenv("NSLOTS"));
pool = parpool(ncores);

addpath('model');

%% model setup
opt = [];
modelClass = [];
rsoa = 100; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.stimContrasts = [.64; .64];
opt.aAI = 0;
opt.aAV = 0;

opt.sigma1 = 0.1;

opt.tauE1 = 100;
opt.tauS1 = 50;

opt.dt = 2;
opt.T = 2.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.stimSet = '10deg';
opt.decoderType = 'continuous';

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

%% generate model output for each parameter combination
seqList = 1:10; % 1 = 0° (identical), 2-10 increasing in 10° steps
suppList = [Inf .04 .1 .2 .4 1 0]; % p, see Eqn 11
contrList = [.64 .64 0; .64 0 .64];

paramList = combvec(suppList,seqList,contrList);

r1_noniden = nan(length(paramList),12,opt.nt);

parfor ii=1:length(paramList)
    opt2 = opt;
    opt2.p_supp = paramList(1,ii);
    opt2.stimContrasts = paramList(3:4,ii);
    iseq = paramList(2,ii);
    [~,p,~] = runModel(opt2,modelClass,400,iseq,rcond);
    r1_noniden(ii,:,:) = p.r1;
end

save('output/respAdapt_noniden.mat','r1_noniden');
