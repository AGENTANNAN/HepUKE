### Paper 2103.01480v2: e+e- -> eta psi(2S) at 4.236-4.600 GeV, 14 energy points
### ORDINARY multi-energy scan
### psi(2S) -> pi+ pi- J/psi, J/psi -> l+ l- (e/mu), eta -> gamma gamma
### 4C kinematic fit chi2 < 40

### 14 energy points from the dataset table
data_4237 = DatasetManager.real_data.find("703_4237")    # 4.236 GeV
data_4245 = DatasetManager.real_data.find("703_4245")    # 4.242 GeV
data_4246 = DatasetManager.real_data.find("703_4246")    # 4.244 GeV
data_4260 = DatasetManager.real_data.find("703_4260")    # 4.258 GeV
data_4270 = DatasetManager.real_data.find("703_4270")    # 4.267 GeV
data_4280 = DatasetManager.real_data.find("703_4280")    # 4.278 GeV
data_4310 = DatasetManager.real_data.find("703_4310")    # 4.308 GeV
data_4360 = DatasetManager.real_data.find("703_4360")    # 4.358 GeV
data_4390 = DatasetManager.real_data.find("703_4390")    # 4.387 GeV
data_4420 = DatasetManager.real_data.find("703_4420")    # 4.416 GeV
data_4470 = DatasetManager.real_data.find("703_4470")    # 4.467 GeV
data_4530 = DatasetManager.real_data.find("703_4530")    # 4.527 GeV
data_4575 = DatasetManager.real_data.find("703_4575")    # 4.575 GeV
data_4600 = DatasetManager.real_data.find("703_4600")    # 4.600 GeV

all_datasets = [data_4237, data_4245, data_4246, data_4260, data_4270,
                data_4280, data_4310, data_4360, data_4390, data_4420,
                data_4470, data_4530, data_4575, data_4600]

### Inclusive MC for the primary energy (4.258 GeV as the reference)
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

### Decay card: e+e- -> eta psi(2S), psi(2S) -> pi+ pi- J/psi, J/psi -> l+ l-
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1  eta  psi(2S)    PHSP;
    Enddecay

    Decay psi(2S)
    1  pi+  pi-  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1  e+  e-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

### Exclusive MC for all 14 energy points (create_exclusive_mc_for)
exMCs = DatasetManager.create_exclusive_mc_for(all_datasets) do |c|
  c.sample_name = "exmc_etapsi2S"
  c.events = 100000
  c.decay_card = decay_card
  c.cross_section = :default
end

### ======================
### Event Selection
### ======================

alg = Algorithm.new("etapsi2S")
alg.set_header(["etapsi2SAlg/etapsi2S.h"])
   .set_constant({"ECMS" => [:double, 4.260]})   # placeholder; multi-energy scan skips ECMS validation

sel = Selection.new
sel.select_track do
     cos_theta 0.93
     Vz 10.0
     Vr 1.0
     nTot "==4"
     nChrp "==2"
     nChrn "==2"
   end
   .select_photon do
     energyThreshold 0.025
     angle_to_track 10.0
     nGam ">=2"    # eta -> gamma gamma
   end
   ### PID: high-momentum leptons (p > 1.0 GeV/c) + pions (p < 0.8 GeV/c)
   .pid(method: :probability) do
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 1.0
     identify :pion, against: [:kaon]
     npip "==1"
     npim "==1"
     nlp "==1"
     nlm "==1"
   end
   ### 4C kinematic fit: gamma gamma pi+ pi- l+ l- with chi2 < 40
   .kinematic_fit([:gamma, :gamma, :pip, :pim, :lp, :lm]) do
     constrain_four_momentum
     chi2_cut 40
     nominal
   end

### Capture inexpressible BOSS-side procedures
alg.note(:jpsi_mass_window, "J/psi mass window: 3064.6 < M(l+l-) < 3140.8 MeV/c^2 applied in ROOT stage")
alg.note(:etap_veto, "eta' -> pi+ pi- eta veto: M(pi+ pi- gamma gamma) > 1.0 GeV/c^2 applied in ROOT stage")
alg.note(:eta_mass_window, "eta mass window: 507.1 < M(gamma gamma) < 579.1 MeV/c^2 applied in ROOT stage")
alg.note(:psi2S_mass_window, "psi(2S) mass window: 3680.3 < M(pi+ pi- J/psi) < 3692.5 MeV/c^2 applied in ROOT stage")
alg.note(:best_combination, "If more than 2 photons, the gamma gamma pi+ pi- l+ l- combination with smallest chi2_4C is retained; selected in ROOT")

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on(all_datasets + [incMC_4260] + exMCs)