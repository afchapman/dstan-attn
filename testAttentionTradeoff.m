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
opt.AVWeights = [1 0];
opt.distributeVoluntary = 0;

soas = 100:100:800;
weights = 0:0.25:1;

model_params = combvec(weights,weights,soas);

perf = nan(2,length(model_params));
parfor ii=1:length(model_params)
    opt2 = opt;
    opt2.AVWeights = model_params(1:2,ii);
    [~,p,~] = runModel(opt2,modelClass,model_params(3,ii),rseq,rcond);
    perf(:,ii) = p.ev;
end

model_out = reshape(perf,2,length(weights),length(weights),length(soas));

figure(1)
for ii=1:2, for jj=1:length(soas)
    subplot(2,length(soas),(ii-1)*length(soas)+jj)
    imagesc(weights,weights,squeeze(model_out(ii,:,:,jj)))
    title(sprintf("T%d, %d ms",ii,soas(jj)))
    clim([min(model_out(ii,:,:,:),[],'all') max(model_out(ii,:,:,:),[],'all')])
    axis square
end, end

figure(2)
for jj=1:length(soas)
    subplot(2,5,jj)
    imagesc(weights,weights,squeeze(sum([.5;.5].*model_out(:,:,:,jj))))
    title(sprintf("Neutral, %d ms",soas(jj)))
    %clim([min(sum([.5;.5].*model_out),[],'all') max(sum([.5;.5].*model_out),[],'all')])
    axis square
end

figure(3)
for jj=1:length(soas)
    subplot(2,5,jj)
    imagesc(weights,weights,squeeze(sum([.75;.25].*model_out(:,:,:,jj))))
    title(sprintf("Cue T1, %d ms",soas(jj)))
%     clim([min(sum([.75;.25].*model_out),[],'all') max(sum([.75;.25].*model_out),[],'all')])
    axis square
end

figure(4)
for jj=1:length(soas)
    subplot(2,5,jj)
    imagesc(weights,weights,squeeze(sum([.25;.75].*model_out(:,:,:,jj))))
    title(sprintf("Cue T2, %d ms",soas(jj)))
%     clim([min(sum([.25;.75].*model_out),[],'all') max(sum([.25;.75].*model_out),[],'all')])
    axis square
end

% T1 fixed at max, function of SOA and T2 gain
figure(5), subplot(121)
plot(weights,squeeze(sum([.75;.25].*model_out(:,end,:,:))))
title('Cue T1, with T1 at max gain')
xlabel('T2 Weight'), ylabel("Expected d'")

% T2 fixed at max, function of SOA and T1 gain
figure(5), subplot(122)
plot(weights,squeeze(sum([.25;.75].*model_out(:,:,end,:))))
title('Cue T2, with T2 at max gain')
xlabel('T1 Weight'), ylabel("Expected d'")

%% add noise to smooth SOA functions
R = 2e3;
model_noise = model_out + normrnd(0,1e-1,[size(model_out) R]);

model_noise_cueT1 = squeeze(sum([.75;.25].*reshape(model_noise,2,[],length(soas),R)));
model_noise_cueT2 = squeeze(sum([.25;.75].*reshape(model_noise,2,[],length(soas),R)));
model_noise_cueN = squeeze(sum([.5;.5].*reshape(model_noise,2,[],length(soas),R)));

[~,idx_cueT1] = max(model_noise_cueT1);
[~,idx_cueT2] = max(model_noise_cueT2);
[~,idx_cueN]  = max(model_noise_cueN);

t1_w = reshape(meshgrid(1:length(weights))',1,[]);
t2_w = reshape(meshgrid(1:length(weights)),1,[]);

max_t1_cueT1 = squeeze(t1_w(idx_cueT1));
max_t2_cueT1 = squeeze(t2_w(idx_cueT1));
max_t1_cueT2 = squeeze(t1_w(idx_cueT2));
max_t2_cueT2 = squeeze(t2_w(idx_cueT2));
max_t1_cueN  = squeeze(t1_w(idx_cueN));
max_t2_cueN  = squeeze(t2_w(idx_cueN));

attn_out = nan(2,3,length(soas));
for ii=1:R, for tt=1:length(soas)
    attn_out(:,1,tt,ii) = model_out(:,max_t1_cueT1(tt,ii),max_t2_cueT1(tt,ii),tt);
    attn_out(:,2,tt,ii) = model_out(:,max_t1_cueT2(tt,ii),max_t2_cueT2(tt,ii),tt);
    attn_out(:,3,tt,ii) = model_out(:,max_t1_cueN(tt,ii), max_t2_cueN(tt,ii), tt);
end, end

attn_outM = mean(attn_out,4);

figure
subplot(121)
plot(soas, squeeze(attn_outM(1,:,:)))
title("T1")

subplot(122)
plot(soas, squeeze(attn_outM(2,:,:)))
title("T2")
