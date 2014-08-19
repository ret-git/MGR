% This script computes causal support for the contingencies in two experiments 
% by Buehner and Cheng (1997)

% generative:

master_obs = [8 8; 6 6; 4 4; 2 2; 0 0;  8 6; 6 4; 4 2; 2 0; 8 ...
              4; 6 2; 4 0; 8 2; 6 0; 8 0];
master_n = [ 8 8; 8 8; 8 8; 8 8; 8 8; 8 8; 8 8; 8 8; 8 ...
             8; 8 8; 8 8; 8 8; 8 8; 8 8; 8 8 ];
n_trials = length(master_n); 

support_gen = zeros(n_trials,1);
for trial = 1:n_trials
  support_gen(trial) = supportsampler(master_obs(trial,:),master_n(trial,:),500000);
end

disp('Generative causal relations:');
disp(' '); 
disp('    A         B         C         D         Support (generative)'); 
disp([master_obs(:,1) master_n(:,1)-master_obs(:,1) master_obs(:,2) master_n(:,2)-master_obs(:,2) support_gen]); 
disp('Hit any key to continue')
pause

% preventive:

master_obs = [8 8; 6 6; 4 4; 2 2; 0 0;  6 8; 4 6; 2 4; 0 2; 4 8; 2 6; 0 4; ...
	     2 8; 0 6; 0 8];
master_n = [ 8 8; 8 8; 8 8; 8 8; 8 8; 8 8; 8 8; 8 8; 8 8; 8 8; ...
	     8 8; 8 8; 8 8; 8 8; 8 8 ];
n_trials = length(master_n); 

support_prev = zeros(n_trials,1);
for trial = 1:n_trials
  support_prev(trial) = supportsampler(master_obs(trial,:),master_n(trial,:),500000,-1);
end

disp(' '); 
disp('Preventive causal relations:');
disp(' '); 
disp('    A         B         C         D         Support (preventive)'); 
disp([master_obs(:,1) master_n(:,1)-master_obs(:,1) master_obs(:,2) master_n(:,2)-master_obs(:,2) support_prev]); 
