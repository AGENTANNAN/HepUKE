# BESIII analysis: search for e+e- -> eta eta psi(2S) and
#                             e+e- -> eta psi_0(4360), psi_0(4360) -> eta psi(2S)
# at Ecm = 4.84, 4.92, 4.95 GeV (0.9 fb^-1 total).
# Partial-reconstruction technique: one eta -> gamma gamma reconstructed,
# the other eta treated as missing. psi(2S) -> pi+ pi- J/psi, J/psi -> e+e-/mu+mu-.

### Datasets ###
data_4840 = DatasetManager.real_data.find("707_4840")   # sqrt(s)=4843.07 MeV
data_4914 = DatasetManager.real_data.find("707_4914")   # sqrt(s)=4918.02 MeV
data_4946 = DatasetManager.real_data.find("707_4946")   # sqrt(s)=4950.93 MeV

incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

all_data  = [data_4840, data_4914, data_4946]
all_incMC = [incMC_4840, incMC_4914, incMC_4946]

### Decay cards ###
# Signal I: e+e- -> eta eta psi(2S)   (both eta's inclusive; only one reconstructed)
# Signal II: e+e- -> eta psi_0(4360), psi_0(4360) -> eta psi(2S)
#   Use psi(4260) as the KKMC top mother (BESIII convention for continuum-like ISR).
# We use two lepton channels for J/psi: ee-mode and mu-mu-mode.
decay_card_etaetapsip_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0000  eta  eta  psi(2S)                   PHSP;
  Enddecay

  Decay psi(2S)
  1.0000  pi+  pi-  J/psi                     JPIPI;
  Enddecay

  Decay J/psi
  1.0000  e+  e-                              VLL;
  Enddecay

  Decay eta
  1.0000  gamma  gamma                        PHSP;
  Enddecay

  End
DECAYCARD

decay_card_etaetapsip_mm = <<~DECAYCARD
  Decay psi(4260)
  1.0000  eta  eta  psi(2S)                   PHSP;
  Enddecay

  Decay psi(2S)
  1.0000  pi+  pi-  J/psi                     JPIPI;
  Enddecay

  Decay J/psi
  1.0000  mu+  mu-                            VLL;
  Enddecay

  Decay eta
  1.0000  gamma  gamma                        PHSP;
  Enddecay

  End
DECAYCARD

# psi_0(4360) is not a standard EvtGen particle: use a psi(2S) alias as placeholder.
decay_card_etapsi04360_ee = <<~DECAYCARD
  Alias psi0_4360  psi(2S)

  Decay psi(4260)
  1.0000  eta  psi0_4360                      PHSP;
  Enddecay

  Decay psi0_4360
  1.0000  eta  psi(2S)                        PHSP;
  Enddecay

  Decay psi(2S)
  1.0000  pi+  pi-  J/psi                     JPIPI;
  Enddecay

  Decay J/psi
  1.0000  e+  e-                              VLL;
  Enddecay

  Decay eta
  1.0000  gamma  gamma                        PHSP;
  Enddecay

  End
DECAYCARD

decay_card_etapsi04360_mm = <<~DECAYCARD
  Alias psi0_4360  psi(2S)

  Decay psi(4260)
  1.0000  eta  psi0_4360                      PHSP;
  Enddecay

  Decay psi0_4360
  1.0000  eta  psi(2S)                        PHSP;
  Enddecay

  Decay psi(2S)
  1.0000  pi+  pi-  J/psi                     JPIPI;
  Enddecay

  Decay J/psi
  1.0000  mu+  mu-                            VLL;
  Enddecay

  Decay eta
  1.0000  gamma  gamma                        PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC ###
exMC_etaetapsip_ee = DatasetManager.create_exclusive_mc_for(all_data) do |c|
  c.sample_name    = "etaeta_psip_ee"
  c.events         = 100_000
  c.decay_card     = decay_card_etaetapsip_ee
  c.cross_section  = :default
end
exMC_etaetapsip_mm = DatasetManager.create_exclusive_mc_for(all_data) do |c|
  c.sample_name    = "etaeta_psip_mm"
  c.events         = 100_000
  c.decay_card     = decay_card_etaetapsip_mm
  c.cross_section  = :default
end
exMC_etapsi04360_ee = DatasetManager.create_exclusive_mc_for([data_4914, data_4946]) do |c|
  c.sample_name    = "eta_psi04360_ee"
  c.events         = 100_000
  c.decay_card     = decay_card_etapsi04360_ee
  c.cross_section  = :default
end
exMC_etapsi04360_mm = DatasetManager.create_exclusive_mc_for([data_4914, data_4946]) do |c|
  c.sample_name    = "eta_psi04360_mm"
  c.events         = 100_000
  c.decay_card     = decay_card_etapsi04360_mm
  c.cross_section  = :default
end

exMC_etaetapsip_ee.each  { |m| m.save_to_config(format: :yaml, file_path: 'exMC_etaetapsip_ee') }
exMC_etaetapsip_mm.each  { |m| m.save_to_config(format: :yaml, file_path: 'exMC_etaetapsip_mm') }
exMC_etapsi04360_ee.each { |m| m.save_to_config(format: :yaml, file_path: 'exMC_etapsi04360_ee') }
exMC_etapsi04360_mm.each { |m| m.save_to_config(format: :yaml, file_path: 'exMC_etapsi04360_mm') }

########################################################################
# Common event selection
#   4 charged tracks, net charge 0.
#   Tracks with p > 1.0 GeV/c => leptons from J/psi; else => pions.
#   Leptons with E_EMC > 1.0 GeV => electrons; E_EMC < 0.4 GeV => muons.
#   At least 2 photons.
#   pi0/eta -> gamma gamma with 1C Kalman fit; require M(gg) in [499.5, 576.9] MeV/c^2.
#   Then partial reconstruction with the second eta as the missing particle.
########################################################################

########################################################################
# Algorithm I : e+e- -> eta eta psi(2S), J/psi -> e+ e-   (ee-mode)
########################################################################
alg_ee = Algorithm.new("EtaEtaPsipEE")
alg_ee.set_header(["EtaEtaPsipEEAlg/EtaEtaPsipEE.h"])
      .set_constant({ "ECMS" => [:double, 4.918] })

sel_ee = Selection.new
sel_ee.select_track {
         cos_theta 0.93           # |cos(theta)| < 0.93
         Vz        10.0           # |Vz| < 10 cm
         Vr        1.0            # |Vxy| < 1 cm
         nChrp     "==2"          # exactly 2 positive
         nChrn     "==2"          # exactly 2 negative
         nNet      "==0"          # net charge = 0
       }
       .select_photon {
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         tdc_emc_start     0
         tdc_emc_end       14
         angle_to_track    10.0
         nGam              ">=2"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]
         identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                        treat_as_electron_if_energy_above: 0.6
         npip ">=1"
         npim ">=1"
       }
       # Kalman fit: reconstruct eta_1 from gamma gamma
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
         chi2_cut 200
         neta ">=1"
       }
       # Nominal kinematic fit: total 4-momentum conservation for
       # e+e- -> eta_1 (reconstructed) + eta_missing + pi+ pi- e+ e- (from psi(2S))
       # with eta_1 mass-constrained (from the Kalman step above)
       .kinematic_fit([:eta, :pip, :pim, :ep, :em]) {
         nominal
         constrain_four_momentum
         miss_track_of(:eta)         # the other eta is missing
         chi2_cut 200
       }

alg_ee
  .note(:lepton_pion_separation,
        "Charged tracks with p > 1.0 GeV/c are assigned as leptons from J/psi; " \
        "others as pions from psi(2S). Leptons with E_EMC > 1.0 GeV are " \
        "identified as electrons; those with E_EMC < 0.4 GeV as muons.")
  .note(:eta_gg_mass_window,
        "Before the 1C Kalman fit, M(gamma gamma) required in " \
        "[499.5, 576.9] MeV/c^2 (Punzi FOM optimised).")
  .note(:eta2_recoil_mass,
        "The second eta is treated as missing; its recoil mass against " \
        "eta_1 psi(2S) is corrected as M(eta_2) = M_recoil + M(l+l-) - M(J/psi). " \
        "Signal region [516.3, 579.5] MeV/c^2.")
  .note(:jpsi_mass_window,
        "M(l+l-) required in [3040.0, 3140.8] MeV/c^2.")
  .note(:psi2S_mass_window,
        "psi(2S) signal region: [3678.3, 3693.5] MeV/c^2 (+/-3 sigma). " \
        "Applied on the corrected M(pi+ pi- J/psi) = M(pi+pi-l+l-) + M(l+l-) - M(J/psi).")
  .with_decay_card(decay_card_etaetapsip_ee)
  .apply(sel_ee)

alg_ee.execute_on(all_data + all_incMC + exMC_etaetapsip_ee)

########################################################################
# Algorithm II : e+e- -> eta eta psi(2S), J/psi -> mu+ mu-   (mumu-mode)
########################################################################
alg_mm = Algorithm.new("EtaEtaPsipMuMu")
alg_mm.set_header(["EtaEtaPsipMuMuAlg/EtaEtaPsipMuMu.h"])
      .set_constant({ "ECMS" => [:double, 4.918] })

sel_mm = Selection.new
sel_mm.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     "==2"
         nChrn     "==2"
         nNet      "==0"
       }
       .select_photon {
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         tdc_emc_start     0
         tdc_emc_end       14
         angle_to_track    10.0
         nGam              ">=2"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]
         identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                        treat_as_electron_if_energy_above: 0.6
         npip ">=1"
         npim ">=1"
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
         chi2_cut 200
         neta ">=1"
       }
       .kinematic_fit([:eta, :pip, :pim, :mup, :mum]) {
         nominal
         constrain_four_momentum
         miss_track_of(:eta)
         chi2_cut 200
       }

alg_mm
  .note(:lepton_pion_separation,
        "Charged tracks with p > 1.0 GeV/c are assigned as leptons from J/psi; " \
        "others as pions from psi(2S). Leptons with E_EMC < 0.4 GeV are " \
        "identified as muons.")
  .note(:muc_depth_cut,
        "For the mu-mu mode, at least one muon candidate is required to have " \
        "hit depth in the MUC greater than 30 cm to suppress pi->mu misID.")
  .note(:eta_gg_mass_window,
        "Before 1C Kalman fit, M(gamma gamma) required in [499.5, 576.9] MeV/c^2.")
  .note(:eta2_recoil_mass,
        "Missing eta recoil-mass signal window [516.3, 579.5] MeV/c^2.")
  .note(:jpsi_mass_window,
        "M(mu+ mu-) required in [3040.0, 3140.8] MeV/c^2.")
  .note(:psi2S_mass_window,
        "psi(2S) signal region [3678.3, 3693.5] MeV/c^2 (+/-3 sigma).")
  .note(:psi04360_mass_window,
        "For the psi_0(4360) search: M(eta_h psi(2S)) window " \
        "[4310.8, 4416.4] MeV/c^2 (+/-3 sigma around 4366 MeV/c^2), where " \
        "eta_h denotes the higher-momentum eta candidate between eta_1 and eta_2.")
  .with_decay_card(decay_card_etaetapsip_mm)
  .apply(sel_mm)

alg_mm.execute_on(all_data + all_incMC + exMC_etaetapsip_mm +
                  exMC_etapsi04360_ee + exMC_etapsi04360_mm)
