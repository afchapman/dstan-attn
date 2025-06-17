%% D-STAN results
%  Example code demonstrating each of the main results from our manuscript
%  Some analyses assume you have downloaded the corresponding .mat files
%  from the associated OSF repository and stored them in an 'output'
%  subfolder.
%  OSF link: https://osf.io/qy9pa/

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

%% reverse correlation
load('output/revCorFits.mat');

% 10 x 10 taus (exc/supp, 0-900 ms) x 601 timepoints (3000 ms at 5 ms spacing)
rc_d1 = reshape(rc_out.d1,10,10,601);
rc_s1 = reshape(rc_out.s1,10,10,601);
rc_r1 = reshape(rc_out.r1,10,10,601);

% figures for tauE = 400, tauS = 100
t = -3000:5:0;

figure
subplot(131)
plot(t,squeeze(rc_d1(5,2,:)))
xlim([-1200 0])

subplot(132)
plot(t,squeeze(rc_s1(5,2,:)))
xlim([-1200 0])

subplot(133)
plot(t,squeeze(rc_r1(5,2,:)))
xlim([-1200 0])

% difference of Gamma fits for all tau combinations
diffOfGamma = @(p,t) p(4,:)'.*(t.*exp(t./p(1,:)') - p(3,:)'.*t.*exp(t./p(2,:)'));

rc_fit = diffOfGamma(rc_out.fitParams,-3000:5:0);
rc_fit = reshape(rc_fit,10,10,601);

figure
subplot(121)
plot(-3000:5:0,squeeze(rc_fit(2:end,5,:)))
xlim([-1200 0])
title('tauE = 400 ms')

subplot(122)
plot(-3000:5:0,squeeze(rc_fit(5,2:end,:)))
xlim([-1200 0])
title('tauS = 400 ms')

%% surround suppression
load('output/surroundSupp.mat');

contrC = [0 logspace(log10(0.05),log10(1),9)]; % center
contrS = [0 .12 .25 .5 1]; % surround

contrResp = reshape(sum(r1_surr(:,6,:),3),length(contrC),length(contrS));

figure
semilogx(contrC,contrResp)
xticks([.05 .1:.1:1])
xlabel('Center contrast')
ylabel('Response')
legend({'Surround: 0%','Surround 12%','Surround 25%','Surround 50%','Surround 100%'},'Location','NorthWest')

%% transient / sustained dynamics

opt.stimContrasts = [.64; 0];
opt.stimDur = 2000;
rsoa = 2001; % avoiding overlap with first stimulus

opt.dt = 2;
opt.T = 8.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

tau_list = [50 100:100:800];
r1_E = nan(12,opt.nt,length(tau_list));
r1_S = nan(12,opt.nt,length(tau_list));

% can parfor this for faster execution if needed - takes a while because of
% the longer analysis window
for ii=1:length(tau_list)
    % varying excitatory tau
    opt2 = opt;
    opt2.tauE1 = tau_list(ii);
    opt2.tauS1 = 0;
    [~,p,~] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    r1_E(:,:,ii) = p.r1;

    % varying suppressive tau
    opt2.tauS1 = tau_list(ii);
    opt2.tauE1 = 0;
    [~,p,~] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    r1_S(:,:,ii) = p.r1;
end

% plot response timecourses
figure
subplot(121)
plot(opt.tlist,squeeze(r1_E(6,:,:)./max(r1_E(6,:,:),[],2)))
axis([0 8500 0 1.2])

subplot(122)
plot(opt.tlist,squeeze(r1_S(6,:,:)./max(r1_S(6,:,:),[],2)))
axis([0 8500 0 1.2])

% time to peak / half max / stable
resp_ttp = nan(2,9);
resp_tthm = nan(1,9);
for ii=1:length(tau_list)
    resp_ttp(1,ii) = opt.tlist(find(r1_E(6,:,ii)>max(r1_E(6,:,ii))*.99,1))-500;
    resp_ttp(2,ii) = opt.tlist(find(r1_S(6,:,ii)==max(r1_S(6,:,ii)),1))-500;

    resp_tthm(ii) = opt.tlist(find(r1_E(6,1250:end,ii)<max(r1_E(6,:,ii))*.5,1));
end

% time to peak
figure
subplot(141)
plot(tau_list,resp_ttp(1,:))
axis([0 850 200 1200])
title('Exc time to peak')

subplot(142)
plot(tau_list,resp_ttp(2,:))
axis([0 850 30 60])
title('Supp time to peak')

% excitatory half max
subplot(143)
plot(tau_list,resp_tthm)
axis([0 850 0 3000])
title('Exc half max')

% suppressive stable level
subplot(144)
plot(tau_list,squeeze(r1_S(6,1250,:)./max(r1_S(6,:,:),[],2)))
axis([0 850 0 0.6])
title('Supp stable level')

%% subadditivity

opt.tauE1 = 100;
opt.tauS1 = 50;

opt.dt = 2;
opt.T = 4.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

% loop over different stim durations and extract sensory response
stimDur = [30 60 120 240 480];
r1_dur = nan(12,opt.nt,length(stimDur));
for ii=1:length(stimDur)
    opt.stimDur = stimDur(ii);
    [~,p,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_dur(:,:,ii) = p.r1;
end

% responses in best neuron
figure
plot(p.tlist,squeeze(r1_dur(6,:,:)))
xlim([0 2000])

% summed response over time for each duration
figure
plot(stimDur,squeeze(sum(r1_dur(6,:,:),2)))
xticks(stimDur)
ylim([0 300])

% effect of taus on subadditivity
load('output/respSubadd.mat');
r1_sub = reshape(r1_sub(:,6,:),10,10,5,3501);

% we have 10 taus [0,50,100:100:800] and 5 stim durations
% divide each response by the previous duration
ra_sub = squeeze(sum(r1_sub(:,:,2:end,:),4)./sum(r1_sub(:,:,1:end-1,:),4));

figure
subplot(121)
plot(stimDur(2:end),squeeze(ra_sub(:,1,:)))
ylim([1 2])
title('tauE')

subplot(122)
plot(stimDur(2:end),squeeze(ra_sub(1,:,:)))
ylim([1 2])
title('tauS')

% plot different parameter combinations
figure
plot(stimDur(2:end),reshape(ra_sub([3 6],[2 6],:),4,4)')
ylim([1 2])
legend({"E = 100, S = 50","E = 100, S = 500","E = 500, S = 50","E = 500, S = 500"})

%% response adaptation
%  first we show response adaptation for repeating an identical stimulus

opt.stimDur = 300;

opt.tauE1 = 100;
opt.tauS1 = 50;

opt.dt = 2;
opt.T = 4.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

soas = [400; 900]; % 100, 600 ISIs
r1_iden = nan(opt.nt,2,2);
for ii=1:length(soas)
    % present both stimuli
    opt.stimContrasts = [.64; .64];

    % rseq=2 presents an identical stimulus
    [~,p,~] = runModel(opt,modelClass,soas(ii),2,rcond);
    r1_iden(:,1,ii) = p.r1(6,:);

    % present only one stimulus
    opt.stimContrasts = [.64; 0];
    [~,p,~] = runModel(opt,modelClass,soas(ii),2,rcond);
    r1_iden(:,2,ii) = p.r1(6,:);
end

% response adaptation in T2 at short SOAs
figure
plot(opt.tlist,r1_iden(:,:,1))
xlim([0 2500])
legend({'T2 present','T2 absent'})

% but not long SOAs
figure
plot(opt.tlist,r1_iden(:,:,2))
xlim([0 2500])
legend({'T2 present','T2 absent'})

% effect of taus on response adaptation
load('output/respAdapt_iden.mat');
r1_iden = reshape(r1_iden(:,6,:),10,10,16,3501);

% we have 10 taus [0,50,100:100:800] and 15 SOAs (+1 for T1 absent)
% see Eqn 10
ra_iden = 1 - sum(r1_iden(:,:,2:end,:)-r1_iden(:,:,1,:),4)./sum(r1_iden(:,:,1,:),4);

figure
subplot(121)
plot(100:100:1500,squeeze(ra_iden(:,1,:)))
ylim([-.1 1])
title('tauE')

subplot(122)
plot(100:100:1500,squeeze(ra_iden(1,:,:)))
ylim([-.1 1])
title('tauS')

% plot different parameter combinations
figure
plot(100:100:1500,reshape(ra_iden([3 6],[2 6],:),4,15))
ylim([-.1 1])
legend({"E = 100, S = 50","E = 100, S = 500","E = 500, S = 50","E = 500, S = 500"})

%% response adaptation in non-identical stimuli
%  testing response adaptation while varying the similarity between stimuli

opt.stimDur = 300;

opt.tauE1 = 100;
opt.tauS1 = 50;

opt.dt = 2;
opt.T = 4.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.stimSet = '10deg';
opt.decoderType = 'continuous';

r1_orth = nan(2,opt.nt,2); % need responses from 2 neurons

% present both stimuli
opt.stimContrasts = [.64; .64];
[~,p,~] = runModel(opt,modelClass,400,10,rcond);
r1_orth(:,:,1) = p.r1([6 12],:);

% present only one stimulus
opt.stimContrasts = [0; .64];
[~,p,~] = runModel(opt,modelClass,400,10,rcond);
r1_orth(:,:,2) = p.r1([6 12],:);

% response adaptation in T2 with orthogonal stimuli
figure
plot(opt.tlist,r1_orth(:,:,1),...
     opt.tlist,r1_orth(:,:,2))
xlim([0 2500])
legend({'T1, T1 present','T2, T1 present','T1, T1 absent','T2, T1 absent'})


% effect of stimulus similarity and suppressive tuning on adaptation
load('output/respAdapt_noniden.mat');

seqList = 1:10; % 1 = 0° (identical), 2-10 increasing in 10° steps
suppList = [Inf .04 .1 .2 .4 1 0]; % p, see Eqn 11
r1_noniden = reshape(r1_noniden(:,12,:),length(suppList),length(seqList),3,[]);

ra_noniden = 1-sum(r1_noniden(:,:,1,:)-r1_noniden(:,:,2,:),4)./sum(r1_noniden(:,:,3,:),4); % Eqn 12

figure
plot(0:10:90,ra_noniden)
ylim([-.1 1])
xlabel('Distance between stimuli (°)')
ylabel('Adaptation index')

%% backward masking

opt.stimDur = 30;

opt.tauE1 = 100;
opt.tauS1 = 50;

soas = [250; 500];
r1_mask = nan(2,opt.nt,2,2);
for ii=1:length(soas)
    % present both stimuli
    opt.stimContrasts = [.64; .64];

    [~,p,~] = runModel(opt,modelClass,soas(ii),3,rcond);
    r1_mask(:,:,1,ii) = p.r1([6 12],:);

    % present only one stimulus
    opt.stimContrasts = [.64; 0];
    [~,p,~] = runModel(opt,modelClass,soas(ii),3,rcond);
    r1_mask(:,:,2,ii) = p.r1([6 12],:);
end

% backward masking in T1 at short SOAs
figure
plot(opt.tlist,r1_mask(:,:,1,1),...
     opt.tlist,r1_mask(:,:,2,1))
xlim([0 2500])
legend({'T1, T1 present','T2, T1 present','T1, T1 absent','T2, T1 absent'})

% but not long SOAs
figure
plot(opt.tlist,r1_mask(:,:,1,2),...
     opt.tlist,r1_mask(:,:,2,2))
xlim([0 2500])
legend({'T1, T1 present','T2, T1 present','T1, T1 absent','T2, T1 absent'})

% effect of taus on backward masking
load('output/backwardMask.mat');
r1_mask = reshape(r1_mask(:,[6 12],:),10,10,16,2,3501);

% we have 10 taus [0,50,100:100:800] and 15 SOAs (+1 for T1 absent)
% see Eqn 13
ra_mask = 1 - sum(r1_mask(:,:,2:end,1,:),5)./sum(r1_mask(:,:,1,2,:),5);

figure
subplot(121)
plot(100:100:1500,squeeze(ra_mask(:,1,:)))
ylim([-.1 1])
title('tauE')

subplot(122)
plot(100:100:1500,squeeze(ra_mask(1,:,:)))
ylim([-.1 1])
title('tauS')

% plot different parameter combinations
figure
plot(100:100:1500,reshape(ra_mask([3 6],[2 6],:),4,15))
ylim([-.1 1])
legend({"E = 100, S = 50","E = 100, S = 500","E = 500, S = 50","E = 500, S = 500"})

%% Contrast-dependent suppression

opt.stimDur = 30;
rsoa = 250;

opt.dt = 2;
opt.T = 4.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.tauE1 = 400;
opt.tauS1 = 100;

r1_supp = nan(2,opt.nt,4);
d_supp = nan(2,4);

% T1 high, T2 high
opt.stimContrasts = [.64; .64];
[~,p,~] = runModel(opt,modelClass,rsoa,3,rcond);
r1_supp(:,:,1) = p.r1([6 12],:);
d_supp(:,1) = p.ev;

% T1 high, T2 low
opt.stimContrasts = [.64; .16];
[~,p,~] = runModel(opt,modelClass,rsoa,3,rcond);
r1_supp(:,:,2) = p.r1([6 12],:);
d_supp(:,2) = p.ev;

% T1 low, T2 high
opt.stimContrasts = [.16; .64];
[~,p,~] = runModel(opt,modelClass,rsoa,3,rcond);
r1_supp(:,:,3) = p.r1([6 12],:);
d_supp(:,3) = p.ev;

% T1 low, T2 low
opt.stimContrasts = [.16; .16];
[~,p,~] = runModel(opt,modelClass,rsoa,3,rcond);
r1_supp(:,:,4) = p.r1([6 12],:);
d_supp(:,4) = p.ev;


% sensory response as a function of T2 contrast
figure
plot(opt.tlist,r1_supp(1,:,1),'r-',...
     opt.tlist,r1_supp(2,:,1),'b-',...
     opt.tlist,r1_supp(1,:,2),'r--',...
     opt.tlist,r1_supp(2,:,2),'b--')
xlim([0 2500])
legend({'T1, T2 high','T2, T2 high','T1, T2 low','T2, T2 low'})

% sensory response as a function of T1 contrast
figure
plot(opt.tlist,r1_supp(1,:,1),'r-',...
     opt.tlist,r1_supp(2,:,1),'b-',...
     opt.tlist,r1_supp(1,:,3),'r--',...
     opt.tlist,r1_supp(2,:,3),'b--')
xlim([0 2500])
legend({'T1, T1 high','T2, T1 high','T1, T1 low','T2, T1 low'})

% model d' as a function of each contrast level
figure
subplot(121), title('T1')
plot(1:2,d_supp(1,[3 1]),'r-',...
     1:2,d_supp(1,[4 2]),'r--')
axis([[.8 2.2] [0 2.5]])
xticks(1:2)
xticklabels({'T1 low','T1 high'})
legend({'T2 high','T2 low'})

subplot(122), title('T2')
plot(1:2,d_supp(2,[2 1]),'b-',...
     1:2,d_supp(2,[4 3]),'b--')
axis([[.8 2.2] [0 2.5]])
xticks(1:2)
xticklabels({'T2 low','T2 high'})
legend({'T1 high','T1 low'})

% contrast-dependent suppression as a function of taus
load('output/contrDepSupp.mat');

figure
subplot(131)
imagesc(0:50:1000,0:50:1000,perf_T1')
axis square
title('T1 modulation')

subplot(132)
imagesc(0:50:1000,0:50:1000,perf_T2')
axis square
title('T2 modulation')

subplot(133)
imagesc(0:50:1000,0:50:1000,perf_T1'.*perf_T2')
axis square
title('Joint modulation')
