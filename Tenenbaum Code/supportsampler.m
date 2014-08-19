function logscore = supportsampler(master_obs, master_n,loops,valence)

% SUPPORTSAMPLER    Computes causal support for contingency data
%
% Usage:
%
%   LOGSCORE = SUPPORTSAMPLER(MASTER_OBS, MASTER_N, LOOPS, VALENCE) 
%
% LOGSCORE    is the log likelihood ratio of a causal relationship
% MASTER_OBS  gives the frequency with which the effect occurs in the
%             presence and absence of the cause (cells A and C)
% MASTER_N    gives the frequency with which the cause was present and
%             absent (A+B and C+D)
% LOOPS       is the number of samples used in the Monte Carlo
%             approximation (default 1000)
% VALENCE     is 1 for generative causes, -1 for preventive causes 
%             (i.e. noisy-OR or noisy-AND-NOT) (default 1)
%
% e.g. for a contingency table [ A B; C D ], use
%           logscore = supportsampler([A C],[A+B C+D],100000,1)

conts = [master_obs(1) master_obs(2) ...
	 master_n(1)-master_obs(1) master_n(2)-master_obs(2)];

if (nargin < 3)
  loops = 1000;
end

if (nargin < 4)
  valence = 1;
end

% assuming a uniform prior on w_C:
power = rand(loops,2);
% for a more general beta prior, use: (requires stats toolbox)
% power = betarnd(ones(loops,2),ones(loops,2));
if (valence > 0)
  probs = [ 1-(1-power(:,1)).*(1-power(:,2)) ...
	    power(:,2) ...
	    (1-power(:,1)).*(1-power(:,2)) ...
	    1-power(:,2)];
else
  probs = [ (1-power(:,1)).*power(:,2) ...
	     power(:,2) ...
	     1-(1-power(:,1)).*power(:,2) ...
	     1-power(:,2)];
end
loglike = sum((ones(loops,1)*conts).*log(probs),2);

logmax = max(loglike);
logscore = logmax + log(sum(exp(loglike-logmax))) - log(loops);

% assuming a uniform prior on w_B
logscore = logscore - betaln(sum(master_obs)+1, ...
			     sum(master_n)-sum(master_obs)+1);

