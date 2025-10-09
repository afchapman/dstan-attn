function response = rfSpatialResponse(theta, x, ntheta, nx, mtheta, mx)

if nargin < 3
    ntheta = 12;
end

if nargin < 4
    nx = 2*length(x)-1;
end

if nargin < 5
    mtheta = 2*ntheta-1;
end

if nargin < 6
    mx = 4*nx-1;
end

responseTh = abs(cos((1:ntheta)*pi/ntheta - theta').^mtheta);

x2 = linspace(2*x(1),2*x(end),nx);
responseX  = abs(cos((x2-x')/range(x2)*2*pi^2/nx).^mx);

% combine responses
response = reshape(responseTh,[],1,ntheta) .* shiftdim(reshape(responseX,[],1,nx),-1);
response = reshape(response,[],ntheta*nx);
