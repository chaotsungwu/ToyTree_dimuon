# A simple dimuon decay reconstruction example on the high energy physics (Monte Carlo Method)
### *You can understand the simple reconstruction strategy from the following steps to execute my codes. All the code are programmed by C and must be executed by CERN ROOT.*
### 1. Use TTree to create a new data file from the following parameter settings (randomly-generated data, maybe 1000000 events):
> double Z_mass = 90 + gRandom->Gaus(0, 5);
> TLorentzVector z(10, 0.0, 0.0, Z_mass); 
> Double_t masses[2] = { 0.106, 0.106};
> TGenPhaseSpace event;
> event.SetDecay(z, 2, masses);
### 2. Mass reconstruction: toytree.C -> toytree_fit.C
