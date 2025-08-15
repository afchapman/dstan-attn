%% SCC parpool

% set up parallel processing
ncores = str2num(getenv('NSLOTS'));
maxNumCompThreads(ncores);

% check if parpool is running
pool = gcp('nocreate');
if isempty(pool)
	myCluster = parcluster('local');
	myCluster.JobStorageLocation = getenv('TMPDIR');
	pool = parpool(myCluster,ncores); % buffer
end

%% model setup
addpath('model');

opt = [];
modelClass = [];
rsoa = []; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 1; % cueT1, cueT2

opt.scaling1 = 5e3;
opt.scaling2 = 5e3;

opt.stimContrasts = [.64; .64];
opt.aAI = 0;
opt.aAV = 2e2;
opt.sigma1 = 0.1;

opt.dt = 2;
opt.T = 3.0*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

opt.tau1 = 52;
opt.tau2 = 52;

opt.tauE1 = 50;
opt.tauS1 = 0;

opt.tauE2 = 400;
opt.tauS2 = 100;

opt.p_supp = Inf;
opt.AVWeights = [1 1];
opt.distributeVoluntary = 0;

%% loop over params of interest
soas = [100:50:500 800];
weights = 0:0.1:1;
taus = [0 50 100:100:800];

model_params = combvec(weights,weights,taus,taus,soas);
perf = nan(2,length(model_params));

% set up parfor loop tracker
q = parallel.pool.DataQueue;
afterEach(q,@parforTracker);
parforTracker(length(model_params),[],1000);

parfor ii=1:length(model_params)
    opt2 = opt;
    opt2.AVWeights = model_params(1:2,ii);
    opt2.tauE2 = model_params(3,ii);
    opt2.tauS2 = model_params(4,ii);
    [~,p,~] = runModel(opt2,modelClass,model_params(5,ii),rseq,rcond);
    perf(:,ii) = p.ev;

    send(q,[]);
end

model_out = reshape(perf,2,length(weights),length(weights),length(taus),length(taus),length(soas));

save('dstan_attn_param_search.mat','model_out');
