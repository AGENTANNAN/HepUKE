# Core DSL classes and dependencies are loaded automatically at execution.
# Analysis: measurement of the tau mass from an energy scan near the tau+tau-
# threshold.  Thirteen two-prong final states are studied:
#   photon-less : ee, e-mu, e-pi, e-K, mu-mu, mu-pi, mu-K, pi-pi, pi-K, K-K
#   X-rho       : e-rho, mu-rho, pi-rho   (rho -> pi pi0, pi0 -> gamma gamma)
# The undetected tau-decay neutrinos are left unreconstructed.

# ============================== Datasets =================================
data_3542 = DatasetManager.real_data.find("710_3542")   # 3542.4 MeV, below threshold -> pure-background control
data_3554 = DatasetManager.real_data.find("710_3554")   # 3553.8 MeV
data_3561 = DatasetManager.real_data.find("710_3561")   # 3561.1 MeV
data_3600 = DatasetManager.real_data.find("710_3600")   # 3600.2 MeV

incMC_3542 = DatasetManager.inclusive_mc.find("710_3542")
incMC_3554 = DatasetManager.inclusive_mc.find("710_3554")
incMC_3561 = DatasetManager.inclusive_mc.find("710_3561")
incMC_3600 = DatasetManager.inclusive_mc.find("710_3600")

scan_points = [data_3542, data_3554, data_3561, data_3600]
scan_incMC  = [incMC_3542, incMC_3554, incMC_3561, incMC_3600]

# ======================= Signal decay cards (EvtGen) =====================
# e+e- -> tau+ tau- ; top mother psi(4260) (KKMC convention).  tau decays use
# TAULNUNU (leptonic) and TAUHADNU (one-prong hadronic); PHOTOS adds FSR.

# --- ten photon-less two-prong channels ---
dc_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 e+ nu_e anti-nu_tau PHOTOS TAULNUNU;
  Enddecay
  Decay tau-
  1.0000 e- anti-nu_e nu_tau PHOTOS TAULNUNU;
  Enddecay
  End
DECAYCARD

dc_emu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 e+ nu_e anti-nu_tau PHOTOS TAULNUNU;
  Enddecay
  Decay tau-
  1.0000 mu- anti-nu_mu nu_tau PHOTOS TAULNUNU;
  Enddecay
  End
DECAYCARD

dc_epi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 e+ nu_e anti-nu_tau PHOTOS TAULNUNU;
  Enddecay
  Decay tau-
  1.0000 pi- nu_tau PHOTOS TAUHADNU;
  Enddecay
  End
DECAYCARD

dc_eK = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 e+ nu_e anti-nu_tau PHOTOS TAULNUNU;
  Enddecay
  Decay tau-
  1.0000 K- nu_tau PHOTOS TAUHADNU;
  Enddecay
  End
DECAYCARD

dc_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 mu+ nu_mu anti-nu_tau PHOTOS TAULNUNU;
  Enddecay
  Decay tau-
  1.0000 mu- anti-nu_mu nu_tau PHOTOS TAULNUNU;
  Enddecay
  End
DECAYCARD

dc_mupi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 mu+ nu_mu anti-nu_tau PHOTOS TAULNUNU;
  Enddecay
  Decay tau-
  1.0000 pi- nu_tau PHOTOS TAUHADNU;
  Enddecay
  End
DECAYCARD

dc_muK = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 mu+ nu_mu anti-nu_tau PHOTOS TAULNUNU;
  Enddecay
  Decay tau-
  1.0000 K- nu_tau PHOTOS TAUHADNU;
  Enddecay
  End
DECAYCARD

dc_pipi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 pi+ nu_tau PHOTOS TAUHADNU;
  Enddecay
  Decay tau-
  1.0000 pi- anti-nu_tau PHOTOS TAUHADNU;
  Enddecay
  End
DECAYCARD

dc_piK = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 pi+ nu_tau PHOTOS TAUHADNU;
  Enddecay
  Decay tau-
  1.0000 K- anti-nu_tau PHOTOS TAUHADNU;
  Enddecay
  End
DECAYCARD

dc_KK = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 K+ nu_tau PHOTOS TAUHADNU;
  Enddecay
  Decay tau-
  1.0000 K- anti-nu_tau PHOTOS TAUHADNU;
  Enddecay
  End
DECAYCARD

# --- three X-rho channels (rho -> pi pi0, pi0 -> gamma gamma) ---
dc_erho = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 e+ nu_e anti-nu_tau PHOTOS TAULNUNU;
  Enddecay
  Decay tau-
  1.0000 rho- nu_tau PHOTOS TAUHADNU;
  Enddecay
  Decay rho-
  1.0000 pi- pi0 VSS;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_murho = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 mu+ nu_mu anti-nu_tau PHOTOS TAULNUNU;
  Enddecay
  Decay tau-
  1.0000 rho- nu_tau PHOTOS TAUHADNU;
  Enddecay
  Decay rho-
  1.0000 pi- pi0 VSS;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_pirho = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  Decay tau+
  1.0000 pi+ nu_tau PHOTOS TAUHADNU;
  Enddecay
  Decay tau-
  1.0000 rho- nu_tau PHOTOS TAUHADNU;
  Enddecay
  Decay rho-
  1.0000 pi- pi0 VSS;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# --- generic tau-pair sample (inclusive tau decays handled by the generator) ---
dc_taupair = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau- PHSP;
  Enddecay
  End
DECAYCARD

# ======================= Exclusive MC (every scan point) =================
exMC_ee = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_ee"
  config.events        = 200_000
  config.decay_card    = dc_ee
  config.cross_section = :default
end

exMC_emu = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_emu"
  config.events        = 200_000
  config.decay_card    = dc_emu
  config.cross_section = :default
end

exMC_epi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_epi"
  config.events        = 200_000
  config.decay_card    = dc_epi
  config.cross_section = :default
end

exMC_eK = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_eK"
  config.events        = 200_000
  config.decay_card    = dc_eK
  config.cross_section = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_mumu"
  config.events        = 200_000
  config.decay_card    = dc_mumu
  config.cross_section = :default
end

exMC_mupi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_mupi"
  config.events        = 200_000
  config.decay_card    = dc_mupi
  config.cross_section = :default
end

exMC_muK = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_muK"
  config.events        = 200_000
  config.decay_card    = dc_muK
  config.cross_section = :default
end

exMC_pipi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_pipi"
  config.events        = 200_000
  config.decay_card    = dc_pipi
  config.cross_section = :default
end

exMC_piK = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_piK"
  config.events        = 200_000
  config.decay_card    = dc_piK
  config.cross_section = :default
end

exMC_KK = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_KK"
  config.events        = 200_000
  config.decay_card    = dc_KK
  config.cross_section = :default
end

exMC_erho = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_erho"
  config.events        = 200_000
  config.decay_card    = dc_erho
  config.cross_section = :default
end

exMC_murho = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_murho"
  config.events        = 200_000
  config.decay_card    = dc_murho
  config.cross_section = :default
end

exMC_pirho = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_pirho"
  config.events        = 200_000
  config.decay_card    = dc_pirho
  config.cross_section = :default
end

exMC_taupair = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_scan_taupair"
  config.events        = 500_000
  config.decay_card    = dc_taupair
  config.cross_section = :default
end

# ==================== Common selection building blocks ===================
# Photon-less two-prong channels: exactly one + and one - track, no photon.
sel_2prong_base = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr         1.0
    nChrp    "==1"
    nChrn    "==1"
    nNet     "==0"
  }
  .select_photon {
    nGam     "==0"      # the ten photon-less channels allow no photon
  }

# X-rho channels: exactly one + and one - track, exactly two good photons.
sel_Xrho_base = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr         1.0
    nChrp    "==1"
    nChrn    "==1"
    nNet     "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end   15        # EMC timing 0-750 ns
    angle_to_track 20.0     # > 20 deg to the nearest charged track
    energyThreshold_b 0.025 # EMC barrel threshold 25 MeV
    energyThreshold_e 0.050 # EMC endcap threshold 50 MeV
    nGam     "==2"
  }

# ================== Photon-less two-prong channels =======================
# ee channel
alg_ee = Algorithm.new("TauMassScanEE")
alg_ee.set_header(["TauMassScanEEAlg/TauMassScanEE.h"])
      .set_constant({"ECMS" => [:double, 3.6002]})
sel_ee = sel_2prong_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  }
  .kinematic_fit([:lp, :lm]) {
    nominal
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_ee.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.60-1.90 GeV for the two-prong channels")
alg_ee.note(:lepton_flavour_separation, "ee / e-mu / mu-mu share the same BOSS lepton PID (index_lp/index_lm combine e and mu); the e vs mu channel separation is performed offline from the NTuple")
alg_ee.with_decay_card(dc_ee).apply(sel_ee)
alg_ee.execute_on(scan_points + scan_incMC + exMC_ee)

# e-mu channel
alg_emu = Algorithm.new("TauMassScanEMu")
alg_emu.set_header(["TauMassScanEMuAlg/TauMassScanEMu.h"])
       .set_constant({"ECMS" => [:double, 3.6002]})
sel_emu = sel_2prong_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  }
  .kinematic_fit([:lp, :lm]) {
    nominal
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_emu.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.60-1.90 GeV for the two-prong channels")
alg_emu.note(:lepton_flavour_separation, "e vs mu separation of the two lepton tracks is performed offline from the NTuple (E/p, MUC depth)")
alg_emu.with_decay_card(dc_emu).apply(sel_emu)
alg_emu.execute_on(scan_points + scan_incMC + exMC_emu)

# e-pi channel
alg_epi = Algorithm.new("TauMassScanEPi")
alg_epi.set_header(["TauMassScanEPiAlg/TauMassScanEPi.h"])
       .set_constant({"ECMS" => [:double, 3.6002]})
sel_epi = sel_2prong_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    nlp "==1"
    npim "==1"
  }
  .kinematic_fit([:lp, :pim]) {
    nominal
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_epi.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.60-1.90 GeV for the two-prong channels")
alg_epi.note(:lepton_flavour_separation, "e-pi vs mu-pi separated offline (same BOSS lepton PID for e and mu)")
alg_epi.with_decay_card(dc_epi).apply(sel_epi)
alg_epi.execute_on(scan_points + scan_incMC + exMC_epi)

# e-K channel
alg_eK = Algorithm.new("TauMassScanEK")
alg_eK.set_header(["TauMassScanEKAlg/TauMassScanEK.h"])
      .set_constant({"ECMS" => [:double, 3.6002]})
sel_eK = sel_2prong_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    identify :kaon, against: [:pion]
    nlp "==1"
    nkm "==1"
  }
  .kinematic_fit([:lp, :km]) {
    nominal
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_eK.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.60-1.90 GeV for the two-prong channels")
alg_eK.note(:lepton_flavour_separation, "e-K vs mu-K separated offline (same BOSS lepton PID for e and mu)")
alg_eK.with_decay_card(dc_eK).apply(sel_eK)
alg_eK.execute_on(scan_points + scan_incMC + exMC_eK)

# mu-mu channel
alg_mumu = Algorithm.new("TauMassScanMuMu")
alg_mumu.set_header(["TauMassScanMuMuAlg/TauMassScanMuMu.h"])
        .set_constant({"ECMS" => [:double, 3.6002]})
sel_mumu = sel_2prong_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  }
  .kinematic_fit([:lp, :lm]) {
    nominal
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_mumu.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.60-1.90 GeV for the two-prong channels")
alg_mumu.note(:lepton_flavour_separation, "mu-mu share the same lepton PID as the other leptonic channels; channel separation done offline from the NTuple")
alg_mumu.with_decay_card(dc_mumu).apply(sel_mumu)
alg_mumu.execute_on(scan_points + scan_incMC + exMC_mumu)

# mu-pi channel
alg_mupi = Algorithm.new("TauMassScanMuPi")
alg_mupi.set_header(["TauMassScanMuPiAlg/TauMassScanMuPi.h"])
        .set_constant({"ECMS" => [:double, 3.6002]})
sel_mupi = sel_2prong_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    nlp "==1"
    npim "==1"
  }
  .kinematic_fit([:lp, :pim]) {
    nominal
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_mupi.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.60-1.90 GeV for the two-prong channels")
alg_mupi.note(:lepton_flavour_separation, "e-pi vs mu-pi separated offline (same BOSS lepton PID for e and mu)")
alg_mupi.with_decay_card(dc_mupi).apply(sel_mupi)
alg_mupi.execute_on(scan_points + scan_incMC + exMC_mupi)

# mu-K channel
alg_muK = Algorithm.new("TauMassScanMuK")
alg_muK.set_header(["TauMassScanMuKAlg/TauMassScanMuK.h"])
       .set_constant({"ECMS" => [:double, 3.6002]})
sel_muK = sel_2prong_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    identify :kaon, against: [:pion]
    nlp "==1"
    nkm "==1"
  }
  .kinematic_fit([:lp, :km]) {
    nominal
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_muK.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.60-1.90 GeV for the two-prong channels")
alg_muK.note(:lepton_flavour_separation, "e-K vs mu-K separated offline (same BOSS lepton PID for e and mu)")
alg_muK.with_decay_card(dc_muK).apply(sel_muK)
alg_muK.execute_on(scan_points + scan_incMC + exMC_muK)

# pi-pi channel
alg_pipi = Algorithm.new("TauMassScanPiPi")
alg_pipi.set_header(["TauMassScanPiPiAlg/TauMassScanPiPi.h"])
        .set_constant({"ECMS" => [:double, 3.6002]})
sel_pipi = sel_2prong_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:pip, :pim]) {
    nominal
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_pipi.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.60-1.90 GeV for the two-prong channels")
alg_pipi.with_decay_card(dc_pipi).apply(sel_pipi)
alg_pipi.execute_on(scan_points + scan_incMC + exMC_pipi)

# pi-K channel
alg_piK = Algorithm.new("TauMassScanPiK")
alg_piK.set_header(["TauMassScanPiKAlg/TauMassScanPiK.h"])
       .set_constant({"ECMS" => [:double, 3.6002]})
sel_piK = sel_2prong_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    identify :kaon, against: [:pion]
    npip "==1"
    nkm "==1"
  }
  .kinematic_fit([:pip, :km]) {
    nominal
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_piK.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.60-1.90 GeV for the two-prong channels")
alg_piK.with_decay_card(dc_piK).apply(sel_piK)
alg_piK.execute_on(scan_points + scan_incMC + exMC_piK)

# K-K channel
alg_KK = Algorithm.new("TauMassScanKK")
alg_KK.set_header(["TauMassScanKKAlg/TauMassScanKK.h"])
      .set_constant({"ECMS" => [:double, 3.6002]})
sel_KK = sel_2prong_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    identify :kaon, against: [:pion]
    nkp "==1"
    nkm "==1"
  }
  .kinematic_fit([:kp, :km]) {
    nominal
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_KK.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.60-1.90 GeV for the two-prong channels")
alg_KK.with_decay_card(dc_KK).apply(sel_KK)
alg_KK.execute_on(scan_points + scan_incMC + exMC_KK)

# ============================ X-rho channels =============================
# e-rho channel
alg_erho = Algorithm.new("TauMassScanERho")
alg_erho.set_header(["TauMassScanERhoAlg/TauMassScanERho.h"])
        .set_constant({"ECMS" => [:double, 3.6002]})
sel_erho = sel_Xrho_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    nlp "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # 1C fit to pi0 mass
    invariant_mass_of(:gamma, :gamma).between(0.1128, 0.1464)              # gamma-gamma mass window
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:lp, :pim, :pi0]) {
    nominal
    invariant_mass_of(:pim, :pi0).between(0.3765, 1.1955)                  # rho (pi pi0) mass window
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_erho.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.55-1.90 GeV for the X-rho channels")
alg_erho.note(:lepton_flavour_separation, "e-rho vs mu-rho separated offline (same BOSS lepton PID for e and mu)")
alg_erho.with_decay_card(dc_erho).apply(sel_erho)
alg_erho.execute_on(scan_points + scan_incMC + exMC_erho)

# mu-rho channel
alg_murho = Algorithm.new("TauMassScanMuRho")
alg_murho.set_header(["TauMassScanMuRhoAlg/TauMassScanMuRho.h"])
         .set_constant({"ECMS" => [:double, 3.6002]})
sel_murho = sel_Xrho_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    nlp "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).between(0.1128, 0.1464)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:lp, :pim, :pi0]) {
    nominal
    invariant_mass_of(:pim, :pi0).between(0.3765, 1.1955)
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_murho.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.55-1.90 GeV for the X-rho channels")
alg_murho.note(:lepton_flavour_separation, "e-rho vs mu-rho separated offline (same BOSS lepton PID for e and mu)")
alg_murho.with_decay_card(dc_murho).apply(sel_murho)
alg_murho.execute_on(scan_points + scan_incMC + exMC_murho)

# pi-rho channel
alg_pirho = Algorithm.new("TauMassScanPiRho")
alg_pirho.set_header(["TauMassScanPiRhoAlg/TauMassScanPiRho.h"])
         .set_constant({"ECMS" => [:double, 3.6002]})
sel_pirho = sel_Xrho_base.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                    treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).between(0.1128, 0.1464)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pip, :pim, :pi0]) {
    nominal
    invariant_mass_of(:pim, :pi0).between(0.3765, 1.1955)
    miss_track_of(:nu)
    constrain_four_momentum
    chi2_cut 200
  }
alg_pirho.note(:missing_mass_constraint, "undetected tau-decay neutrinos handled by a missing-mass (recoil-mass) window 1.55-1.90 GeV for the X-rho channels")
alg_pirho.with_decay_card(dc_pirho).apply(sel_pirho)
alg_pirho.execute_on(scan_points + scan_incMC + exMC_pirho)

# NOTE: the channel-dependent PTEM and acoplanarity windows are applied offline
# on the NTuple and are therefore outside the BOSS selection scope.