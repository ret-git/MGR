function predictions = causalgenerative(master_obs,master_n,data,bootflag)

% CAUSALGENERATIVE   Computes causal support for a generative relation
%
% Usage:
%
%   PREDICTIONS = CAUSALGENERATIVE(MASTER_OBS,MASTER_N,DATA,BOOTFLAG)
%
% PREDICTIONS is a cellarray struct
% MASTER_OBS  gives the frequency with which the effect occurs in the
%             presence and absence of the cause (cells A and C)
% MASTER_N    gives the frequency with which the cause was present and
%             absent (A+B and C+D)
% DATA        human judgments for contingencies in MASTER_OBS and MASTER_N
% BOOTFLAG    whether to compute a bootstrap confidence interval on the
%             correlation coefficient (default = 0)
%

if (nargin < 4)
  bootflag = 0;
end

n = length(data);

master_p = master_obs./master_n; 
master_deltap = master_p(:,1)-master_p(:,2);
master_cheng = master_deltap./(1-master_p(:,2)+eps);
master_chi2 = zeros(n,1);
master_loglr = zeros(n,1);
master_pci = zeros(n,1);
warning off
for trial = 1:n 
  master_bayes(trial) = supportsampler(master_obs(trial,:),master_n(trial,:),500000);
  D(1,2,2) = master_obs(trial,2);
  D(2,2,2) = master_obs(trial,1);
  D(1,2,1) = master_n(trial,2)-D(1,2,2); 
  D(2,2,1) = master_n(trial,1)-D(2,2,2); 
  D(:,1,:) = 0;
  Nij = squeeze(D(:,2,:))+eps^2; 
  Nidot = sum(Nij,2); 
  Ndotj = sum(Nij,1); 
  N = sum(sum(Nij)); 
  nij = (Nidot*Ndotj)/N;
  chi2 = sum(sum(((Nij-nij).^2)./nij));
  loglr = betaln(Nij(1,1)+1,Nij(1,2)+1)+betaln(Nij(2,1)+1,Nij(2,2)+1)- ...
	  betaln(Ndotj(1)+1,Ndotj(2)+1);
  master_loglr(trial) = loglr;
  master_chi2(trial) = chi2;
  master_pci(trial) = (D(2,2,2)+D(1,2,1)-D(2,2,1)-D(1,2,2)) ... 
      / sum(sum(squeeze(D(:,2,:))));
end
warning on

predictions{1}.name = '\Delta P    ';
predictions{1}.gamma = fminsearch('powertransform',1,[],master_deltap,data);
predictions{1}.value = sign(master_deltap).*abs(master_deltap).^predictions{1}.gamma;
predictions{1}.r = -powertransform(predictions{1}.gamma,master_deltap,data);
predictions{1}.raw = master_deltap;
if bootflag
  predictions{1}.std = bootcheck(predictions{1}.gamma,master_deltap,data);
end

predictions{2}.name = 'Power       ';
predictions{2}.gamma = fminsearch('powertransform',1,[],master_cheng,data);
predictions{2}.value = sign(master_cheng).*abs(master_cheng).^predictions{2}.gamma;
predictions{2}.r = -powertransform(predictions{2}.gamma,master_cheng,data);
predictions{2}.raw = master_cheng;
if bootflag
  predictions{2}.std = bootcheck(predictions{2}.gamma,master_cheng,data);
end

predictions{3}.name = 'pCI         ';
predictions{3}.gamma = fminsearch('powertransform',1,[],master_pci,data);
predictions{3}.value = sign(master_pci).*abs(master_pci).^predictions{3}.gamma;
predictions{3}.r = -powertransform(predictions{3}.gamma,master_pci,data);
predictions{3}.raw = master_pci;
if bootflag
  predictions{3}.std = bootcheck(predictions{3}.gamma,master_pci,data);
end

master_bayes = reshape(master_bayes,length(master_bayes),1);

predictions{4}.name = 'Support (OR)';
predictions{4}.gamma = fminsearch('powertransform',1,[],master_bayes,data);
predictions{4}.value = sign(master_bayes).*abs(master_bayes).^predictions{4}.gamma;
predictions{4}.r = -powertransform(predictions{4}.gamma,master_bayes,data);
predictions{4}.raw = master_bayes;
if bootflag
  predictions{4}.std = bootcheck(predictions{4}.gamma,master_bayes,data);
end

predictions{5}.name = '\chi^2      ';
predictions{5}.gamma = fminsearch('powertransform',1,[],master_chi2,data);
predictions{5}.value = master_chi2.^predictions{5}.gamma;
predictions{5}.r = -powertransform(predictions{5}.gamma,master_chi2,data);
predictions{5}.raw = master_chi2;
if bootflag
  predictions{5}.std = bootcheck(predictions{5}.gamma,master_chi2,data);
end

predictions{6}.name = 'No theory   ';
predictions{6}.gamma = fminsearch('powertransform',1,[],master_loglr, ...
				  data);
predictions{6}.value = sign(master_loglr).*abs(master_loglr).^predictions{6}.gamma;
predictions{6}.r = -powertransform(predictions{6}.gamma,master_loglr,data);
predictions{6}.raw = master_loglr;
if bootflag
  predictions{6}.std = bootcheck(predictions{6}.gamma,master_loglr,data);
end


disp('Generative causal relations:');
disp(' '); 
disp('Model            Fit       Best power transformation');
for i = 1:6
  disp([predictions{i}.name ' :   ' num2str(predictions{i}.r) ...
	'   ' num2str(predictions{i}.gamma,'%1.5f')]);
end
