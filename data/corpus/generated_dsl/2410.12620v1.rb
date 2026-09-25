# =============================================================================
# BOSS DSL — search for e+e- -> phi chi_c0(3415) and phi eta_c2(1D), phi -> K+K-
# at 16 energy points from 4.47 to 4.95 GeV (6.7 fb^-1, BOSS 703/706/707)
# =============================================================================

### Dataset preparation ###
# --- 16 real-data energy points (BOSS 703/706/707) ---
data_703 = %w[703_4470 703_4530 703_4575 703_4600].map { |n| DatasetManager.real_data.find(n) }
data_706 = %w[706_4610 706_4620 706_4640 706_4660 706_4680 706_4700].map { |n| DatasetManager.real_data.find(n) }
data_707 = %w[707_4740 707_4750 707_4780 707_4840 707_4914 707_4946].map { |n| DatasetManager.real_data.find(n) }
data_points = data_703 + data_706 + data_707

# --- corresponding inclusive MC at every point ---
incMC_703 = %w[703_4470 703_4530 703_4575 703_4600].map { |n| DatasetManager.inclusive_mc.find(n) }
incMC_706 = %w[706_4610 706_4620 706_4640 706_4660 706_4680 706_4700].map { |n| DatasetManager.inclusive_mc.find(n) }
incMC_707 = %w[707_4740 707_4750 707_4780 707_4840 707_4914 707_4946].map { |n| DatasetManager.inclusive_mc.find(n) }
incMC_points = incMC_703 + incMC_706 + incMC_707

# --- eta_c2(1D) search restricted to the three points at/above 4.84 GeV ---
eta_c2_data_points  = [DatasetManager.real_data.find("707_4840"),
                       DatasetManager.real_data.find("707_4914"),
                       DatasetManager.real_data.find("707_4946")]
eta_c2_incMC_points = [DatasetManager.inclusive_mc.find("707_4840"),
                       DatasetManager.inclusive_mc.find("707_4914"),
                       DatasetManager.inclusive_mc.find("707_4946")]

### Decay cards (EvtGen format, psi(4260) as top mother for KKMC) ###
# chi_c0 -> pi+ pi-
decay_card_chi_c0_pipi = <<~DECAYCARD
    Decay psi(4260)
    1.0 phi chi_c0 PHSP;
    Enddecay
    Decay phi
    1.0 K+ K- VSS;
    Enddecay
    Decay chi_c0
    1.0 pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# chi_c0 -> pi+ pi- pi0 pi0
decay_card_chi_c0_pipi_pi0pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0 phi chi_c0 PHSP;
    Enddecay
    Decay phi
    1.0 K+ K- VSS;
    Enddecay
    Decay chi_c0
    1.0 pi+ pi- pi0 pi0 PHSP;
    Enddecay
    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# chi_c0 -> K+ K- pi+ pi-
decay_card_chi_c0_kkpipi = <<~DECAYCARD
    Decay psi(4260)
    1.0 phi chi_c0 PHSP;
    Enddecay
    Decay phi
    1.0 K+ K- VSS;
    Enddecay
    Decay chi_c0
    1.0 K+ K- pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# chi_c0 -> 2(pi+ pi-)
decay_card_chi_c0_4pi = <<~DECAYCARD
    Decay psi(4260)
    1.0 phi chi_c0 PHSP;
    Enddecay
    Decay phi
    1.0 K+ K- VSS;
    Enddecay
    Decay chi_c0
    1.0 pi+ pi- pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# chi_c0 -> 3(pi+ pi-)
decay_card_chi_c0_6pi = <<~DECAYCARD
    Decay psi(4260)
    1.0 phi chi_c0 PHSP;
    Enddecay
    Decay phi
    1.0 K+ K- VSS;
    Enddecay
    Decay chi_c0
    1.0 pi+ pi- pi+ pi- pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# eta_c2(1D) -> pi+ pi-
decay_card_eta_c2_pipi = <<~DECAYCARD
    Decay psi(4260)
    1.0 phi eta_c2 PHSP;
    Enddecay
    Decay phi
    1.0 K+ K- VSS;
    Enddecay
    Decay eta_c2
    1.0 pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

### Exclusive MC (50k events per available energy point) ###
# chi_c0 modes are generated at every one of the 16 points
exMC_chi_c0_pipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_phi_chi_c0_pipi"
  config.events        = 50_000
  config.decay_card    = decay_card_chi_c0_pipi
  config.cross_section = :default
end

exMC_chi_c0_pipi_pi0pi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_phi_chi_c0_pipi_pi0pi0"
  config.events        = 50_000
  config.decay_card    = decay_card_chi_c0_pipi_pi0pi0
  config.cross_section = :default
end

exMC_chi_c0_kkpipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_phi_chi_c0_kkpipi"
  config.events        = 50_000
  config.decay_card    = decay_card_chi_c0_kkpipi
  config.cross_section = :default
end

exMC_chi_c0_4pi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_phi_chi_c0_4pi"
  config.events        = 50_000
  config.decay_card    = decay_card_chi_c0_4pi
  config.cross_section = :default
end

exMC_chi_c0_6pi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_phi_chi_c0_6pi"
  config.events        = 50_000
  config.decay_card    = decay_card_chi_c0_6pi
  config.cross_section = :default
end

# eta_c2 mode only at the three points >= 4.84 GeV
exMC_eta_c2_pipi = DatasetManager.create_exclusive_mc_for(eta_c2_data_points) do |config|
  config.sample_name   = "exmc_phi_eta_c2_pipi"
  config.events        = 50_000
  config.decay_card    = decay_card_eta_c2_pipi
  config.cross_section = :default
end

### Event selection (BOSS) ###

# -----------------------------------------------------------------------
# Mode 1: chi_c0 -> pi+ pi-   (final state K+K-pi+pi- = 2K2pi)
# -----------------------------------------------------------------------
alg_chi_c0_pipi = Algorithm.new("PhiChiC0ToPiPi")
alg_chi_c0_pipi.set_header(["PhiChiC0ToPiPiAlg/PhiChiC0ToPiPi.h"])
               .set_constant({"ECMS" => [:double, 4.6]})
               .set_alias({"std::vector<double>" => "Vdouble"})
sel_chi_c0_pipi = Selection.new
  .select_track {              # charged tracks: |cosθ|<0.93, |Vz|<10 cm, Vr<1 cm, 2K2pi
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
    nNet  "==0"
  }
  .select_photon {             # photons: E>25/50 MeV, TDC [0,14], >10 deg from tracks, <4
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    nGam "<4"
  }
  .pid(method: :probability) { # prob PID with 0.001 cut, K/pi separation
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim]) {  # 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 40
  }
alg_chi_c0_pipi.note(:helix_correction, "charged-track helix-parameter correction applied before the 4C kinematic fit; efficiency difference with/without the correction is estimated by re-running the BOSS selection on signal MC")
               .with_decay_card(decay_card_chi_c0_pipi)
               .apply(sel_chi_c0_pipi)
alg_chi_c0_pipi.execute_on(data_points + incMC_points + exMC_chi_c0_pipi)

# -----------------------------------------------------------------------
# Mode 2: chi_c0 -> pi+ pi- pi0 pi0   (K+K-pi+pi- + 2pi0 = 2K2pi)
#   pi0 reconstructed from gamma-gamma (1C Kalman fit, chi2<200) before 4C
# -----------------------------------------------------------------------
alg_chi_c0_pipi_pi0pi0 = Algorithm.new("PhiChiC0ToPiPiPi0Pi0")
alg_chi_c0_pipi_pi0pi0.set_header(["PhiChiC0ToPiPiPi0Pi0Alg/PhiChiC0ToPiPiPi0Pi0.h"])
                       .set_constant({"ECMS" => [:double, 4.6]})
                       .set_alias({"std::vector<double>" => "Vdouble"})
sel_chi_c0_pipi_pi0pi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
    nNet  "==0"
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    nGam "<8"                  # pi0 mode: allow up to 8 photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # constrain each gamma-gamma pair to pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=2"                 # at least two pi0
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :pi0, :pi0]) {  # 4C fit on full final state
    nominal
    constrain_four_momentum
    chi2_cut 40
  }
alg_chi_c0_pipi_pi0pi0.note(:helix_correction, "charged-track helix-parameter correction applied before the 4C kinematic fit; efficiency difference with/without the correction is estimated by re-running the BOSS selection on signal MC")
                       .with_decay_card(decay_card_chi_c0_pipi_pi0pi0)
                       .apply(sel_chi_c0_pipi_pi0pi0)
alg_chi_c0_pipi_pi0pi0.execute_on(data_points + incMC_points + exMC_chi_c0_pipi_pi0pi0)

# -----------------------------------------------------------------------
# Mode 3: chi_c0 -> K+ K- pi+ pi-   (phi K+K- + K+K-pi+pi- = 4K2pi)
# -----------------------------------------------------------------------
alg_chi_c0_kkpipi = Algorithm.new("PhiChiC0ToKKPiPi")
alg_chi_c0_kkpipi.set_header(["PhiChiC0ToKKPiPiAlg/PhiChiC0ToKKPiPi.h"])
                 .set_constant({"ECMS" => [:double, 4.6]})
                 .set_alias({"std::vector<double>" => "Vdouble"})
sel_chi_c0_kkpipi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==3"
    nChrn "==3"
    nNet  "==0"
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    nGam "<4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==2"
    nkm "==2"
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:kp, :kp, :km, :km, :pip, :pim]) {  # 4C fit, 2K+ 2K- pi+ pi-
    nominal
    constrain_four_momentum
    chi2_cut 40
  }
alg_chi_c0_kkpipi.note(:helix_correction, "charged-track helix-parameter correction applied before the 4C kinematic fit; efficiency difference with/without the correction is estimated by re-running the BOSS selection on signal MC")
                 .with_decay_card(decay_card_chi_c0_kkpipi)
                 .apply(sel_chi_c0_kkpipi)
alg_chi_c0_kkpipi.execute_on(data_points + incMC_points + exMC_chi_c0_kkpipi)

# -----------------------------------------------------------------------
# Mode 4: chi_c0 -> 2(pi+ pi-)   (phi K+K- + 2pi+2pi- = 2K4pi)
# -----------------------------------------------------------------------
alg_chi_c0_4pi = Algorithm.new("PhiChiC0To4Pi")
alg_chi_c0_4pi.set_header(["PhiChiC0To4PiAlg/PhiChiC0To4Pi.h"])
              .set_constant({"ECMS" => [:double, 4.6]})
              .set_alias({"std::vector<double>" => "Vdouble"})
sel_chi_c0_4pi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==3"
    nChrn "==3"
    nNet  "==0"
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    nGam "<4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
    identify :pion, against: [:kaon, :proton]
    npip "==2"
    npim "==2"
  }
  .kinematic_fit([:kp, :km, :pip, :pip, :pim, :pim]) {  # 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 40
  }
alg_chi_c0_4pi.note(:helix_correction, "charged-track helix-parameter correction applied before the 4C kinematic fit; efficiency difference with/without the correction is estimated by re-running the BOSS selection on signal MC")
              .with_decay_card(decay_card_chi_c0_4pi)
              .apply(sel_chi_c0_4pi)
alg_chi_c0_4pi.execute_on(data_points + incMC_points + exMC_chi_c0_4pi)

# -----------------------------------------------------------------------
# Mode 5: chi_c0 -> 3(pi+ pi-)   (phi K+K- + 3pi+3pi- = 2K6pi)
# -----------------------------------------------------------------------
alg_chi_c0_6pi = Algorithm.new("PhiChiC0To6Pi")
alg_chi_c0_6pi.set_header(["PhiChiC0To6PiAlg/PhiChiC0To6Pi.h"])
              .set_constant({"ECMS" => [:double, 4.6]})
              .set_alias({"std::vector<double>" => "Vdouble"})
sel_chi_c0_6pi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==4"
    nChrn "==4"
    nNet  "==0"
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    nGam "<4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
    identify :pion, against: [:kaon, :proton]
    npip "==3"
    npim "==3"
  }
  .kinematic_fit([:kp, :km, :pip, :pip, :pip, :pim, :pim, :pim]) {  # 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 40
  }
alg_chi_c0_6pi.note(:helix_correction, "charged-track helix-parameter correction applied before the 4C kinematic fit; efficiency difference with/without the correction is estimated by re-running the BOSS selection on signal MC")
              .with_decay_card(decay_card_chi_c0_6pi)
              .apply(sel_chi_c0_6pi)
alg_chi_c0_6pi.execute_on(data_points + incMC_points + exMC_chi_c0_6pi)

# -----------------------------------------------------------------------
# Mode 6: partial rec — phi chi_c0 with one missing K+-
#   chi_c0 -> pi+ pi-, one K from phi undetected => 1C fit, chi2<45
# -----------------------------------------------------------------------
alg_chi_c0_pipi_missK = Algorithm.new("PhiChiC0ToPiPiMissK")
alg_chi_c0_pipi_missK.set_header(["PhiChiC0ToPiPiMissKAlg/PhiChiC0ToPiPiMissK.h"])
                     .set_constant({"ECMS" => [:double, 4.6]})
                     .set_alias({"std::vector<double>" => "Vdouble"})
sel_chi_c0_pipi_missK = Selection.new
  .select_track {              # reconstructed: K- pi+ pi-  (K+ from phi is missing)
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==1"
    nChrn "==2"
    nNet  "==-1"
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    nGam "<4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkm "==1"
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:pip, :pim, :km]) {     # 1C fit with the missing K+
    nominal
    miss_track_of :kp
    constrain_four_momentum
    chi2_cut 45
  }
alg_chi_c0_pipi_missK.with_decay_card(decay_card_chi_c0_pipi)
                     .apply(sel_chi_c0_pipi_missK)
alg_chi_c0_pipi_missK.execute_on(data_points + incMC_points + exMC_chi_c0_pipi)

# -----------------------------------------------------------------------
# Mode 7: partial rec — chi_c0 -> pi+ pi- pi0 pi0 with one missing pi0
#   reconstruct K+K-pi+pi- + one pi0, miss one pi0 => 1C fit, chi2<6
# -----------------------------------------------------------------------
alg_chi_c0_pipi_pi0pi0_missPi0 = Algorithm.new("PhiChiC0ToPiPiPi0Pi0MissPi0")
alg_chi_c0_pipi_pi0pi0_missPi0.set_header(["PhiChiC0ToPiPiPi0Pi0MissPi0Alg/PhiChiC0ToPiPiPi0Pi0MissPi0.h"])
                              .set_constant({"ECMS" => [:double, 4.6]})
                              .set_alias({"std::vector<double>" => "Vdouble"})
sel_chi_c0_pipi_pi0pi0_missPi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
    nNet  "==0"
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    nGam "<8"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct the single visible pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :pi0]) {  # 1C fit with the missing pi0
    nominal
    miss_track_of :pi0
    constrain_four_momentum
    chi2_cut 6
  }
alg_chi_c0_pipi_pi0pi0_missPi0.with_decay_card(decay_card_chi_c0_pipi_pi0pi0)
                              .apply(sel_chi_c0_pipi_pi0pi0_missPi0)
alg_chi_c0_pipi_pi0pi0_missPi0.execute_on(data_points + incMC_points + exMC_chi_c0_pipi_pi0pi0)

# -----------------------------------------------------------------------
# Mode 8: eta_c2(1D) -> pi+ pi-   (phi K+K- + pi+pi- = 2K2pi)  4C, chi2<60
# -----------------------------------------------------------------------
alg_eta_c2_pipi = Algorithm.new("PhiEtaC2ToPiPi")
alg_eta_c2_pipi.set_header(["PhiEtaC2ToPiPiAlg/PhiEtaC2ToPiPi.h"])
               .set_constant({"ECMS" => [:double, 4.85]})
               .set_alias({"std::vector<double>" => "Vdouble"})
sel_eta_c2_pipi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
    nNet  "==0"
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    nGam "<4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim]) {  # 4C fit, chi2<60
    nominal
    constrain_four_momentum
    chi2_cut 60
  }
alg_eta_c2_pipi.note(:helix_correction, "charged-track helix-parameter correction applied before the 4C kinematic fit; efficiency difference with/without the correction is estimated by re-running the BOSS selection on signal MC")
               .with_decay_card(decay_card_eta_c2_pipi)
               .apply(sel_eta_c2_pipi)
alg_eta_c2_pipi.execute_on(eta_c2_data_points + eta_c2_incMC_points + exMC_eta_c2_pipi)