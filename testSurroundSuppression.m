
addpath('model');

% default model settings
opt = [];
modelClass = [];
rsoa = 2001; % SOA = 250 ms (see runModel)
rseq = 2; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.aAI = 0;
opt.aAV = 0;
opt.sigma1 = 0.1;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

opt.tauE1 = 500;
opt.tauS1 = 500;

opt.dt = 2;
opt.T = 2.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.stimDur = 100;

opt.x = [-1 1];
opt.nx = 2;

opt.stimMode = 'surr_supp';

contrC = [0 exp(linspace(log(0.1),log(1),19))];
contrS = [0 .12 .25 .5 1];

paramList = combvec(contrC,contrS);

all_d1 = nan(24,opt.nt,length(paramList));
all_s1 = nan(24,opt.nt,length(paramList));
all_f1 = nan(24,opt.nt,length(paramList));
all_r1 = nan(24,opt.nt,length(paramList));
for ii=1:length(paramList)
    opt.stimContrasts = paramList(:,ii);
    [~,p,~] = runModel(opt, modelClass, rsoa, rseq, rcond);
    all_d1(:,:,ii) = p.d1;
    all_s1(:,:,ii) = p.s1;
    all_f1(:,:,ii) = p.f1;
    all_r1(:,:,ii) = p.r1;
end

contrResp = reshape(sum(all_r1(6,:,:),2),length(contrC),length(contrS));

figure
plot(contrC,contrResp)
xlabel('Center contrast')
ylabel('Response')
legend({'Surround: 0%','Surround 12%','Surround 25%','Surround 50%','Surround 100%'},'Location','NorthWest')
