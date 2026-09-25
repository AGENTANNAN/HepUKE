### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC

# --- Decay cards: psi(3686) -> pi0 h_c ; h_c -> five modes ---
# Mode I: h_c -> p pbar pi+ pi-
decay_card_pppipi = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.000 p+ anti-p- pi+ pi- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode II: h_c -> pi+ pi- pi0
decay_card_pipipi0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.000 pi+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode III: h_c -> 2(pi+ pi-) pi0
decay_card_2pipipi0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.000 pi+ pi- pi+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode IV: h_c -> 3(pi+ pi-) pi0
decay_card_3pipipi0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.000 pi+ pi- pi+ pi- pi+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode V: h_c -> K+ K- pi+ pi-
decay_card_kkpipi = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.000 K+ K- pi+ pi- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# --- Exclusive MC (100k events each) ---
exMC_pppipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_pppipi"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_pppipi
  config.cross_section   = :default
end

exMC_pipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_pipipi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_pipipi0
  config.cross_section   = :default
end

exMC_2pipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_2pipipi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_2pipipi0
  config.cross_section   = :default
end

exMC_3pipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_3pipipi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_3pipipi0
  config.cross_section   = :default
end

exMC_kkpipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_kkpipi"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_kkpipi
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ===================== Mode I : h_c -> p pbar pi+ pi- =====================
alg_pppipi = Algorithm.new("Pi0HcToPPbarPiPi")
alg_pppipi.set_header(["Pi0HcToPPbarPiPiAlg/Pi0HcToPPbarPiPi.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_pppipi = Selection.new
  .select_track {
    cos_theta 0.93     # |cos(theta)| < 0.93
    Vz        10.0     # |Vz| < 10 cm
    Vr        1.0      # Vr < 1 cm
    nChrp     "==2"    # 2 positive tracks (p+, pi+)
    nChrn     "==2"    # 2 negative tracks (pbar-, pi-)
    nNet      "==0"    # net charge 0
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025   # E > 25 MeV (barrel)
    energyThreshold_e 0.050   # E > 50 MeV (endcap)
    nGam              ">=2"   # at least 2 photons (tag pi0)
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]  # p and pbar
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # remove identified (anti)protons
  .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks are pi+/pi-
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct the tag pi0 from photon pairs
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :prp, :prm, :pip, :pim]) {   # 4C fit to tag pi0 + charged final state
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .kinematic_fit([:gamma, :prp, :prm, :pip, :pim]) { # competing hypothesis psi(3686) -> gamma chi_c2
    constrain_four_momentum
  }
alg_pppipi.with_decay_card(decay_card_pppipi).apply(sel_pppipi)

# ===================== Mode II : h_c -> pi+ pi- pi0 =====================
alg_pipipi0 = Algorithm.new("Pi0HcToPiPiPi0")
alg_pipipi0.set_header(["Pi0HcToPiPiPi0Alg/Pi0HcToPiPiPi0.h"])
           .set_constant({"ECMS" => [:double, 3.686]})
           .note(:tag_pi0_definition,
                 "two pi0's are present; the tag pi0 is taken as the lower-energy pi0 from " \
                 "psi(3686)->pi0 h_c, the other pi0 belongs to the h_c decay")

sel_pipipi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"   # 4 photons -> two pi0's
  }
  .assign({:chrgp => :pip, :chrgn => :pim})    # pure-pion mode: no PID
  .kalman_kinematic_fit([:gamma, :gamma]) {     # reconstruct both pi0's (tag + h_c)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:pi0, :pi0, :pip, :pim]) {    # 4C fit to pi0(pi0) + charged
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .kinematic_fit([:gamma, :pi0, :pip, :pim]) {  # competing hypothesis psi(3686) -> gamma chi_c2
    constrain_four_momentum
  }
alg_pipipi0.with_decay_card(decay_card_pipipi0).apply(sel_pipipi0)

# ===================== Mode III : h_c -> 2(pi+ pi-) pi0 =====================
alg_2pipipi0 = Algorithm.new("Pi0HcTo2PiPiPi0")
alg_2pipipi0.set_header(["Pi0HcTo2PiPiPi0Alg/Pi0HcTo2PiPiPi0.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .note(:tag_pi0_definition,
                  "two pi0's are present; the tag pi0 is the lower-energy pi0 from " \
                  "psi(3686)->pi0 h_c, the other pi0 belongs to the h_c decay")

sel_2pipipi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"
  }
  .assign({:chrgp => :pip, :chrgn => :pim})     # pure-pion mode: no PID
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:pi0, :pi0, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .kinematic_fit([:gamma, :pi0, :pip, :pip, :pim, :pim]) {
    constrain_four_momentum
  }
alg_2pipipi0.with_decay_card(decay_card_2pipipi0).apply(sel_2pipipi0)

# ===================== Mode IV : h_c -> 3(pi+ pi-) pi0 =====================
alg_3pipipi0 = Algorithm.new("Pi0HcTo3PiPiPi0")
alg_3pipipi0.set_header(["Pi0HcTo3PiPiPi0Alg/Pi0HcTo3PiPiPi0.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .note(:tag_pi0_definition,
                  "two pi0's are present; the tag pi0 is the lower-energy pi0 from " \
                  "psi(3686)->pi0 h_c, the other pi0 belongs to the h_c decay")

sel_3pipipi0 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==3"
    nChrn     "==3"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"
  }
  .assign({:chrgp => :pip, :chrgn => :pim})     # pure-pion mode: no PID
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:pi0, :pi0, :pip, :pip, :pip, :pim, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .kinematic_fit([:gamma, :pi0, :pip, :pip, :pip, :pim, :pim, :pim]) {
    constrain_four_momentum
  }
alg_3pipipi0.with_decay_card(decay_card_3pipipi0).apply(sel_3pipipi0)

# ===================== Mode V : h_c -> K+ K- pi+ pi- =====================
alg_kkpipi = Algorithm.new("Pi0HcToKKPiPi")
alg_kkpipi.set_header(["Pi0HcToKKPiPiAlg/Pi0HcToKKPiPi.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_kkpipi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"    # K+, pi+
    nChrn     "==2"    # K-, pi-
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"   # tag pi0
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]  # K+ and K-
    nkp ">=1"
    nkm ">=1"
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])      # remove identified K+/K-
  .assign({:chrgp => :pip, :chrgn => :pim})    # remaining tracks are pi+/pi-
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pi0, :kp, :km, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .kinematic_fit([:gamma, :kp, :km, :pip, :pim]) {   # competing hypothesis psi(3686) -> gamma chi_c2
    constrain_four_momentum
  }
alg_kkpipi.with_decay_card(decay_card_kkpipi).apply(sel_kkpipi)

### Execute ###
root_files_pppipi   = alg_pppipi.execute_on([psip_data, psip_incMC, exMC_pppipi])
root_files_pipipi0  = alg_pipipi0.execute_on([psip_data, psip_incMC, exMC_pipipi0])
root_files_2pipipi0 = alg_2pipipi0.execute_on([psip_data, psip_incMC, exMC_2pipipi0])
root_files_3pipipi0 = alg_3pipipi0.execute_on([psip_data, psip_incMC, exMC_3pipipi0])
root_files_kkpipi   = alg_kkpipi.execute_on([psip_data, psip_incMC, exMC_kkpipi])