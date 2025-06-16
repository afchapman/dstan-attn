
addpath('model');

% default model settings
opt = [];
modelClass = [];
rsoa = 250; % SOA = 250 ms (see runModel) - we'll adjust as needed for each simulation
rseq = []; % default (orthogonal) orientation sequence
rcond = 3; % neutral cues only

opt.aAI = 0; % turn off involuntary attention layer
opt.aAV = 0; % and voluntary attention layer
opt.sigma1 = 0.1; % semi-saturation constant for single sensory layer

opt.scaling1 = 1e4; % scaling for T1
opt.scaling2 = 1e4; % and T2

opt.display.plotTS = 0; % turn off plotting for each simulation
opt.display.plotPerf = 0;

opt.stimDur = 300;

opt.stimSet = '10deg';
opt.decoderType = 'continuous';

opt.tauE1 = 100;
opt.tauS1 = 50;

opt.dt = 2;
opt.T = 2.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.stimContrasts = [.64; .64];

soas = 400; % 100, 600 ISIs

seqList = 1:10;
suppList = [Inf .04 .1 .2 .4 1 0];
paramList = combvec(1:length(seqList),1:length(suppList));

all_r1 = nan(12,opt.nt,length(paramList),2);
base_r1 = nan(12,opt.nt,length(suppList));
for ii=1:length(paramList)
    this_opt = opt;
    iseq = seqList(paramList(1,ii));
    this_opt.m_supp = suppList(paramList(2,ii));

    this_opt.stimContrasts = [.64; .64];
    [~,p,~] = runModel(this_opt,modelClass,soas,iseq,rcond);
    all_r1(:,:,ii,1) = p.r1(12,:);

    this_opt.stimContrasts = [.64; 0];
    [~,p,~] = runModel(this_opt,modelClass,soas,iseq,rcond);
    all_r1(:,:,ii,2) = p.r1(12,:);

    if paramList(1,ii)==1
        this_opt.stimContrasts = [0; .64];
        [~,p,~] = runModel(this_opt,modelClass,soas,iseq,rcond);
        base_r1(:,:,paramList(2,ii)) = p.r1;
    end
end

all_r1 = reshape(all_r1,12,opt.nt,length(seqList),length(suppList),2);

adapt_r1 = squeeze(all_r1(12,:,:,:,1)-all_r1(12,:,:,:,2));

adapt_ind = 1-squeeze(sum(adapt_r1,1) ./ sum(base_r1(12,:,:),2));

figure
plot(0:10:90,adapt_ind)
ylim([-.1 1])
xlabel('Distance between stimuli (°)')
ylabel('Adaptation index')
