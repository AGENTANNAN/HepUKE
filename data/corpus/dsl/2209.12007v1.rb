# Paper: 2209.12007v1
# Title: Search for direct formation of X(3872) in e+e- -> pi+ pi- J/psi
#        via two-photon fusion
# Energy: 3807.7, 3867.4, 3871.3, 3896.2 MeV (4 energy points)
# J/psi -> e+e- and J/psi -> mu+mu- (two independent lepton modes)
# Cross section measurement; upper limit on Gamma_ee x B
# Total ~322 pb^-1

### Dataset preparation ###
# 4 scan energy points around X(3872)
data_709_3686 = DatasetManager.load_real_data.find("709_3686")

all_data = [data_709_3686]
all_incMC = DatasetManager.load_inclusive_mc

# Decay card: e+e- -> pi+ pi- J/psi (continuum, no resonance)
# J/psi -> e+e- for electron mode
decay_card_ee = <<~DECAYCARD
    Decay e+ e-
    1.000  pi+  pi-  J/psi                         PHSP;
    Enddecay

    Decay J/psi
    1.000  e+  e-                                  PHSP;
    Enddecay
End
DECAYCARD

# J/psi -> mu+ mu- for muon mode
decay_card_mumu = <<~DECAYCARD
    Decay e+ e-
    1.000  pi+  pi-  J/psi                         PHSP;
    Enddecay

    Decay J/psi
    1.000  mu+  mu-                                PHSP;
    Enddecay
End
DECAYCARD

exMC_ee = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "ee_to_pipi_jpsi_ee"
  config.events         = 500000
  config.decay_card     = decay_card_ee
  config.cross_section  = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "ee_to_pipi_jpsi_mumu"
  config.events         = 500000
  config.decay_card     = decay_card_mumu
  config.cross_section  = :default
end

### Event selection for J/psi -> e+e- mode ###
alg_ee = Algorithm.new("ee_to_pipi_jpsi_ee")
alg_ee.set_header(["eePiPiJpsiEEAlg/eePiPiJpsiEE.h"])

event_selection_ee = Selection.new

# 4 charged tracks: pi+ pi- e+ e-
event_selection_ee.select_track {
                   cos_theta 0.93
                   Vz   10.0
                   Vr   1.0
                 }
                 .pid(method: :probability) {
                   prob_cut 0.001
                   identify :electron, against: [:pion]
                   identify :pion, against: [:electron, :kaon]
                 }
                 .assign({:pip => :pip, :pim => :pim, :ep => :ep, :em => :em})

# 4C kinematic fit
event_selection_ee.kinematic_fit {
                   constrain_four_momentum
                   nominal
                   chi2_cut 200
                 }

alg_ee.with_decay_card(decay_card_ee).apply(event_selection_ee)

# Event selection details for electron mode
alg_ee.note(:ee_selection,
  "J/psi -> e+e- mode. Lepton tracks: p > 1.0 GeV/c. Pion tracks: p < 0.6 GeV/c. Electron PID: E/p > 1.1 GeV (EMC energy deposit). 4C kinematic fit (chi2<60). Photon conversion veto: cos theta(pi+pi-) < 0.95 and cos theta(pi+- e-+) < 0.98. Signal from fits to M(e+e-) recoil mass of pi+pi- in [3.08,3.12] GeV/c2. Applied in ROOT.")

### Event selection for J/psi -> mu+ mu- mode ###
alg_mumu = Algorithm.new("ee_to_pipi_jpsi_mumu")
alg_mumu.set_header(["eePiPiJpsiMuMuAlg/eePiPiJpsiMuMu.h"])

event_selection_mumu = Selection.new

# 4 charged tracks: pi+ pi- mu+ mu-
event_selection_mumu.select_track {
                     cos_theta 0.93
                     Vz   10.0
                     Vr   1.0
                   }
                   .pid(method: :probability) {
                     prob_cut 0.001
                     identify :muon, against: [:pion]
                     identify :pion, against: [:kaon]
                   }
                   .assign({:pip => :pip, :pim => :pim, :mup => :mup, :mum => :mum})

# 4C kinematic fit
event_selection_mumu.kinematic_fit {
                     constrain_four_momentum
                     nominal
                     chi2_cut 200
                   }

alg_mumu.with_decay_card(decay_card_mumu).apply(event_selection_mumu)

# Event selection details for muon mode
alg_mumu.note(:mumu_selection,
  "J/psi -> mu+mu- mode. Same track momentum cuts as electron mode (p_pi < 0.6, p_mu > 1.0 GeV/c). Muon PID: E/p < 0.35 GeV (EMC energy deposit). Same 4C kinematic fit, photon conversion vetoes. Signal from fits to M(mu+mu-) recoil mass of pi+pi-. Applied in ROOT.")

# Cross section measurement
alg_ee.note(:cross_section,
  "Cross sections measured at 4 c.m. energies: 3807.7 (50.5 pb^-1), 3867.4 (108.9 pb^-1), 3871.3 (110.3 pb^-1), 3896.2 (52.6 pb^-1) MeV. No enhancement at X(3872) peak. Upper limit Gamma_ee x B(X(3872)->pi+pi-J/psi) < 7.5e-3 eV at 90% CL for Gamma_tot=1.19 MeV. Factor ~17 improvement over previous limit. Applied in ROOT.")

all_datasets = all_data + all_incMC + exMC_ee + exMC_mumu
root_files_ee = alg_ee.execute_on(all_datasets)
root_files_mumu = alg_mumu.execute_on(all_datasets)