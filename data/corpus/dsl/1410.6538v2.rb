# BESIII arxiv 1410.6538v2: e+e- -> omega chi_cJ (J=0,1,2), sqrt(s)=4.21-4.42 GeV
# omega -> pi+ pi- pi0; chi_c0 -> pi+pi- / K+K-; chi_c1,2 -> gamma J/psi, J/psi -> l+l-

### Datasets — 9 CME points from 4.21 to 4.42 GeV ###
energies = ["4210", "4220", "4230", "4245", "4260", "4310", "4360", "4390", "4420"]
data_samples   = energies.map { |e| DatasetManager.real_data.find("703_#{e}") }
incMC_samples  = energies.map { |e| DatasetManager.inclusive_mc.find("703_#{e}") }

### Decay cards ###
decay_card_wchic0 = <<~DECAYCARD
  Decay psi(4260)
  0.5000 omega chi_c0                             PHSP;
  0.5000 omega chi_c0                             PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0                               OMEGA_DALITZ;
  Enddecay

  Decay chi_c0
  0.500 pi+ pi-                                   PHSP;
  0.500 K+  K-                                    PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                               PHSP;
  Enddecay

  End
DECAYCARD

decay_card_wchic1 = <<~DECAYCARD
  Decay psi(4260)
  1.000 omega chi_c1                              PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0                               OMEGA_DALITZ;
  Enddecay

  Decay chi_c1
  1.000 gamma J/psi                               PHSP;
  Enddecay

  Decay J/psi
  0.500 e+ e-                                     PHOTOS VLL;
  0.500 mu+ mu-                                   PHOTOS VLL;
  Enddecay

  Decay pi0
  1.000 gamma gamma                               PHSP;
  Enddecay

  End
DECAYCARD

decay_card_wchic2 = <<~DECAYCARD
  Decay psi(4260)
  1.000 omega chi_c2                              PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0                               OMEGA_DALITZ;
  Enddecay

  Decay chi_c2
  1.000 gamma J/psi                               PHSP;
  Enddecay

  Decay J/psi
  0.500 e+ e-                                     PHOTOS VLL;
  0.500 mu+ mu-                                   PHOTOS VLL;
  Enddecay

  Decay pi0
  1.000 gamma gamma                               PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC samples for each mode at each energy ###
exMC_wchic0 = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "omega_chic0"
  config.events        = 200000
  config.decay_card    = decay_card_wchic0
  config.cross_section = :default
end
exMC_wchic1 = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "omega_chic1"
  config.events        = 200000
  config.decay_card    = decay_card_wchic1
  config.cross_section = :default
end
exMC_wchic2 = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "omega_chic2"
  config.events        = 200000
  config.decay_card    = decay_card_wchic2
  config.cross_section = :default
end

### Algorithm 1: omega chi_c0 ###
# Final state: pi+ pi- pi0 (omega) + pi+pi- or K+K- (chi_c0)
alg_wchic0 = Algorithm.new("OmegaChic0")
alg_wchic0.set_header(["OmegaChic0Alg/OmegaChic0.h"])
          .set_constant({"ECMS" => [:double, 4.230]})

sel_wchic0 = Selection.new
sel_wchic0.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     "==2"
             nChrn     "==2"
             nNet      "==0"
          }
          .select_photon {
             nGam              ">=2"
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             tdc_emc_start     0
             tdc_emc_end       14
          }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]
             identify :kaon, against: [:pion, :proton]
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 200
             npi0 ">=1"
          }
          # 5C kinematic fit: total 4-momentum + M(gamma gamma) constrained to m(pi0).
          # Two hypotheses: pi+pi- chi_c0 and K+K- chi_c0; choose the one with smaller chi2_5C.
          .kinematic_fit([:pi0, :pip, :pim, :pip, :pim]) {
             nominal
             constrain_four_momentum
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 100
          }
          .kinematic_fit([:pi0, :pip, :pim, :kp, :km]) {
             constrain_four_momentum
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          }

alg_wchic0.note(:momentum_assignment,
  "Tracks with p > 1 GeV/c are identified as originating from chi_c0; " \
  "lower-momentum pions are assigned to omega decays.")
          .note(:pi_vs_K_hypothesis,
  "For omega chi_c0: two hypotheses (pi+pi- and K+K-) are compared via " \
  "chi2_5C. The lower chi2 hypothesis defines the mode. Stored as separate " \
  "kinematic_fit chi2 values for ROOT-level ranking.")
          .note(:helix_correction,
  "Helix parameter correction applied to charged tracks before the 5C fit.")

alg_wchic0.with_decay_card(decay_card_wchic0).apply(sel_wchic0)
alg_wchic0.execute_on(data_samples + incMC_samples + exMC_wchic0)

### Algorithm 2: omega chi_c1 (and chi_c2), chi_cJ -> gamma J/psi, J/psi -> l+l- ###
# Final state: pi+ pi- pi0 (omega) + gamma + l+ l- (J/psi)
# One extra photon compared with chi_c0.
alg_wchic1 = Algorithm.new("OmegaChic1")
alg_wchic1.set_header(["OmegaChic1Alg/OmegaChic1.h"])
          .set_constant({"ECMS" => [:double, 4.360]})

sel_wchic1 = Selection.new
sel_wchic1.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     "==2"
             nChrn     "==2"
             nNet      "==0"
          }
          .select_photon {
             nGam              ">=3"
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             tdc_emc_start     0
             tdc_emc_end       14
          }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]
             identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                            treat_as_electron_if_energy_above: 0.6
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 200
             npi0 ">=1"
          }
          # 5C kinematic fit: 4-momentum + M(gg)->m(pi0); final state has extra gamma.
          .kinematic_fit([:pi0, :gamma, :pip, :pim, :index_lp, :index_lm]) {
             nominal
             constrain_four_momentum
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 60
          }

alg_wchic1.note(:lepton_ID_by_ECL,
  "Charged particle with ECL deposition > 1 GeV identified as electron; " \
  "otherwise as muon. Handled by identify_high_momentum_leptons.")
          .note(:psi_prime_veto,
  "Events with M(pi+pi- l+l-) or M_recoil(pi+pi-) in [3.68, 3.70] GeV/c^2 " \
  "are vetoed to suppress e+e- -> pi+pi- psi(2S) and pi0pi0 psi(2S) " \
  "backgrounds (applied at ROOT level).")

alg_wchic1.with_decay_card(decay_card_wchic1).apply(sel_wchic1)
alg_wchic1.execute_on(data_samples + incMC_samples + exMC_wchic1)

### Algorithm 3: omega chi_c2 — identical selection to omega chi_c1 (same final state) ###
alg_wchic2 = Algorithm.new("OmegaChic2")
alg_wchic2.set_header(["OmegaChic2Alg/OmegaChic2.h"])
          .set_constant({"ECMS" => [:double, 4.360]})

sel_wchic2 = Selection.new
sel_wchic2.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     "==2"
             nChrn     "==2"
             nNet      "==0"
          }
          .select_photon {
             nGam              ">=3"
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             tdc_emc_start     0
             tdc_emc_end       14
          }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]
             identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                            treat_as_electron_if_energy_above: 0.6
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 200
             npi0 ">=1"
          }
          .kinematic_fit([:pi0, :gamma, :pip, :pim, :index_lp, :index_lm]) {
             nominal
             constrain_four_momentum
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 60
          }

alg_wchic2.note(:psi_prime_veto,
  "Same psi(2S) veto as in omega chi_c1 channel.")

alg_wchic2.with_decay_card(decay_card_wchic2).apply(sel_wchic2)
alg_wchic2.execute_on(data_samples + incMC_samples + exMC_wchic2)
