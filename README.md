# A simple dimuon decay reconstruction example on the high energy physics (Toy Monte Carlo)
### *You can understand the simple reconstruction strategy from the following steps to execute my codes. All the code are programmed by C and must be executed by CERN ROOT.*
### 1. Use TTree to create a new data file called "ToyTree.root" from the following parameter settings (randomly-generated data, maybe 1000000 events)
> double Z_mass = 90 + gRandom -> Gaus(0, 5);\
> TLorentzVector z(10, 0.0, 0.0, Z_mass);\
> Double_t masses[2] = {0.106, 0.106};\
> TGenPhaseSpace event;\
> event.SetDecay(z, 2, masses);
### 2. Mass reconstruction
toytree.C --> toytree_fit.C *(You will derive the mass window of the signal used in the **sideband subtraction** below)*
### 3. pT distribution (2 methods): *"Sideband Subtraction" & "Fitting Method"*
For the **sideband subtraction**: toytree_pt.C --> toytree_final.C\
For the **fitting method**: toytree_final_another.C --> toytree_final_another_fit.C
### 4. Combine the results of the pT distribution from two methods to compare if they are consistent. Also, base on the result of "fitting method" (theoritically, it has better precision than another method), and put the result of "sideband subtraction" on the numerator, you will get the pT distribution of their ratio which is one way to see how far they are different (of course, error analysis is needed to be considered during the calculation).
toytree_final_combine.C --> toytree_final_ratio.C
### *All my result can be found in the "Toy Monte Carlo_Dimuon Decay.pptx", and you can have a check-up after you finish it. Please don't hesitate to let me know if there is any mistakes.*
