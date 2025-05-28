function p = setStim(p)

if strcmp(p.stimMode,'standard')
    % start and end times of T1
    stimStart = p.stimOnset;
    stimEnd = p.stimOnset + p.stimDur;

    % make stim
    timeSeries = zeros([p.norient p.nt]);
    timeSeries(p.stimseq(1),unique(round((stimStart:p.dt:stimEnd)/p.dt))) = 1; % T1
    timeSeries(p.stimseq(2),unique(round(((stimStart:p.dt:stimEnd) + p.soa)/p.dt))) = 1; % T2

    % stimulus contrasts
    contrSeries = zeros([p.norient p.nt]);
    contrSeries(p.stimseq(1),unique(round((stimStart:p.dt:stimEnd)/p.dt))) = p.contrast(1); % T1
    contrSeries(p.stimseq(2),unique(round(((stimStart:p.dt:stimEnd) + p.soa)/p.dt))) = p.contrast(2); % T2

    p.stim = timeSeries;
    p.stimContrast = contrSeries;
elseif strcmp(p.stimMode,'surr_supp') % surround suppression simulation
    stimStart = p.stimOnset;
    stimEnd = p.stimOnset + p.stimDur;

    % make stim
    timeSeries = zeros([p.norient.*p.nx p.nt]);
    timeSeries(p.stimseq(1),unique(round((stimStart:p.dt:stimEnd)/p.dt))) = 1; % center
    timeSeries(p.stimseq(2)+p.norient,unique(round((stimStart:p.dt:stimEnd)/p.dt))) = 1; % surround

    % stimulus contrasts
    contrSeries = zeros([p.norient.*p.nx p.nt]);
    contrSeries(p.stimseq(1),unique(round((stimStart:p.dt:stimEnd)/p.dt))) = p.contrast(1); % center
    contrSeries(p.stimseq(2)+p.norient,unique(round((stimStart:p.dt:stimEnd)/p.dt))) = p.contrast(2); % surround

    p.stim = timeSeries;
    p.stimContrast = contrSeries;
else % random stim array
    p.stim = zeros(4,length(p.tlist));
    p.stim(1,:) = randi([0 1],1,length(p.tlist));
    p.stimContrast = p.stim;
end
