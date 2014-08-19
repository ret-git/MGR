
This folder contains Matlab code for simulating models of human causal
induction, used in Griffiths & Tenenbaum (2005).  

The folder contains eight files:

README			This file

supportsampler.m	Computes causal support via Monte Carlo simulation

support_bc97.m          Script applying supportsampler.m to contingencies from 
                          Buehner and Cheng (1997).  

causalgenerative.m	Compares six models for generative causal relations

causalpreventive.m	Compares six models for preventive causal relations

causal_bc97.m        	Script applying causalgenerative.m and causalpreventive.m
                          to contingencies and behavioral data from B&C (1997)

powertransform.m	Power transformation used in optimizing model fits

bootcheck.m             Bootstrap code for computing confidence intervals on model 
                          fits


The code in causalgenerative.m and causalpreventive.m compares these six models: 

\Delta P 	      	The difference in the probability of the effect
       			in the presence and absence of the cause

Causal Power		As defined in Cheng (1997)

pCI			The proportion of Confirming Instances (White, 2002c)

Causal Support		Bayesian inference with either the Noisy-OR or the
       			Noisy-AND-NOT, depending on valence of causal
       			relation (Griffiths & Tenenbaum, 2005)

\chi^2			The frequentist chi-square test

No theory 		Bayesian inference using the generic
   			parameterization, and uniform priors on the
   			parameters of the model


All references appear in the bibliography of

Griffiths, T. L., & Tenenbaum, J. B. (2005). Structure and strength in
causal induction. Cognitive Psychology.

