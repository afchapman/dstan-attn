function p = setTaskSpatiotemporal(p)

w = p.AVWeights;
attnTemp = p.AVTime;
attnSpace = p.AVSpace;

% check we have matching sizes for temporal and spatial attention
nTemp = size(attnTemp,1);
nSpace = size(attnSpace,1);

if nTemp~=nSpace
    error('temporal and spatial attention have different numbers of stimuli!')
end

% check stimuli match weights
if nTemp==length(w) & nSpace==length(w)
    1;
elseif nTemp~=length(w)
    error("temporal attention times don't match weights");
elseif nSpace~=length(w)
    error("spatial attention positions don't match weights");
end

% make voluntary attention array
timeSeries = zeros(p.ntheta,p.nx,p.nt,nTemp);
x2 = linspace(2*p.x(1),2*p.x(end),p.nx);
for nT=1:nTemp
    times = attnTemp(nT,1):p.dt:attnTemp(nT,2);
    spatFilt = normpdf(x2,p.x(attnSpace(nT,1)),p.AVSpaceWid(nT,1)) .* ...
        sqrt(2*pi) .* p.AVSpaceWid(nT,1) .* w(nT);
    timeSeries(:,:,round(times./p.dt),nT) = repmat(spatFilt,p.ntheta,1,length(times));
end
timeSeries = max(timeSeries,[],4);
timeSeries = reshape(timeSeries,[],p.nt);

p.task = timeSeries;
