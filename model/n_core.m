function [r, f, s] = n_core(d, sigma, p, r_prev, tau, temporalW, s_tuning, dt)

% function [r, f, s] = n_core(d, sigma, p, r_prev, tau, dt)
%
% r_prev = p.r(:,idx-1)
% d = d
% s = p.s(:,idx)
% f = p.f(:,idx)
% r = p.r(:,idx)

% We let the suppressive drive be broad in feature space by letting all
% features contribute equally to the normalization pool

% Normalization pool
pool = abs(d); % abs in case drive is negative

% Suppressive Drive
% s = sum((s_tuning*pool).*temporalW,2);
if any(s_tuning)
    s = sum(s_tuning * (pool.*temporalW),2);
else
    s = sum(pool.*temporalW,'all');
end

% Normalization
f = d(:,end) ./ (s + halfExp(sigma, p));

% Update firing rates
r = r_prev + (dt/tau)*(-r_prev + f);
