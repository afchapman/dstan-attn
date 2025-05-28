function p = n_model(p)

% function p = n_model(p)
% This function is called by runModel.m

idx = 1; % corresponds to t=0
for t = p.dt:p.dt:p.T
    idx = idx+1;
    
    %% Sensory layer 1 (S1)
    % Input (stimulus)
    inp = p.stim(:,idx);
    contrast = p.stimContrast(:,idx);
    
    % Excitatory drive
    if any(inp)
        drive1(:,idx) = halfExp(p.rfresp(logical(inp),:)'*contrast(logical(inp)),p.p); % select pre-calculated response
    else
        drive1(:,idx) = zeros(p.ntheta.*p.nx,1);
    end
    
    attGain = halfExp(1+p.rav(:,idx-1)*p.aAV).*halfExp(1+p.rai(:,idx-1)*p.aAI);
    p.d1(:,idx) = attGain.*sum(drive1.*p.tempWE1(idx:-1:1),2);
    
    % Normalize and update firing rates
    [p.r1(:,idx), p.f1(:,idx), p.s1(:,idx)] = n_core(...
        p.d1(:,1:idx), p.sigma1, p.p, p.r1(:,idx-1), p.tau1, p.tempWS1(idx:-1:1), p.dt);
    
    %% Sensory layer 2 (S2)
    % Excitatory drive
    drive2 = halfExp(p.r1(:,1:idx),p.p);
    p.d2(:,idx) = sum(drive2.*p.tempWE2(idx:-1:1),2);
    
    % Normalize and update firing rates
    [p.r2(:,idx), p.f2(:,idx), p.s2(:,idx)] = n_core(...
        p.d2(:,1:idx), p.sigma2, p.p, p.r2(:,idx-1), p.tau2, p.tempWS2(idx:-1:1), p.dt);
    
    %% Sensory layer 3 (S3)    
    % Excitatory drive
    drive = halfExp(p.r2(:,idx),p.p);
    p.d3(:,idx) = drive;
    
    % Normalize and update firing rates
    [p.r3(:,idx), p.f3(:,idx), p.s3(:,idx)] = n_core(...
        p.d3(:,1:idx), p.sigma3, p.p, p.r3(:,idx-1), p.tau3, 1, p.dt);
    
    %% Decision layer
    % Excitatory drive
    % select response from which to decode
    switch p.modelClass
        case 'LC'
            response = p.r3;
        otherwise
            response = p.r1;
    end
    % decode just between CCW/CW for the appropriate axis
    for iStim = 1:2
        switch p.stimseq(iStim)
            case {1, 2}
                rfresp(:,:,iStim) = p.rfresp(1:2,:);
            case {3, 4}
                rfresp(:,:,iStim) = p.rfresp(3:4,:);
        end
        evidence(iStim) = decodeEvidence(response(1:p.ntheta,idx)', rfresp(:,1:p.ntheta,iStim));
        evidence(iStim) = evidence(iStim)*p.decisionWindows(iStim,idx); % only accumulate if in the decision window
        % evidence(:,abs(evidence)<1e-3) = 0; % otherwise near-zero response will give a little evidence

        % % drive
        % drive = evidence;
        % p.dd(iStim,idx) = drive;
    end

%     evidence(evidence==min(evidence)) = 0; % only accumulate the stronger evidence
    p.dd(:,idx) = evidence;

    % Normalize and update firing rates
    [p.rd(:,idx), p.fd(:,idx), p.sd(:,idx)] = n_core(...
        p.dd(:,idx), p.sigmaD, p.p, p.rd(:,idx-1), p.tauD, 1, p.dt);
    
    %% Voluntary attention layer
    % Inputs
    inp = p.task(:,1:idx-1);
    
    % Excitatory drive
    drive = halfExp(inp, p.p);
    p.dav(:,idx) = sum(drive.*p.tempWEAV(idx-1:-1:1),2); % not feature-specific
    
    % Normalize and update firing rates
    [p.rav(:,idx), p.fav(:,idx), p.sav(:,idx)] = n_core(...
        p.dav(:,1:idx), p.sigmaA, p.p, p.rav(:,idx-1), p.tauAV, p.tempWSAV(idx:-1:1), p.dt);
    
    %% Involuntary attention layer
    % Excitatory drive
    drive = prefilter(p.r1, p.h, p.p, p.dt, idx);
    p.dai(:,idx) = sum(drive); % not feature-specific
    
    % Normalize and update firing rates
    [p.rai(:,idx), p.fai(:,idx), p.sai(:,idx)] = n_core(...
        p.dai(:,idx), p.sigmaA, p.p, p.rai(:,idx-1), p.tauAI, 1, p.dt);
    
end
