# ============================================================================
# e+e- -> pi+ pi- h_c , h_c -> gamma eta_c  at sqrt(s) = 4.23, 4.26, 4.36, 4.42 GeV
# h_c tagged via the recoil mass against pi+pi- in [3.515, 3.535] GeV/c^2
# eta_c reconstructed in four exclusive modes (KKpi0, KS Kpi, 2(pi+pi-pi0), p pbar)
# plus one inclusive (recoil only) mode.  ECMS varies per energy point.
# ============================================================================

### ------------------------- Datasets (four energy points) ------------------------- ###
data_4230  = DatasetManager.real_data.find("703_4230")      # sqrt(s) ~ 4.23 GeV
data_4260  = DatasetManager.real_data.find("703_4260")      # sqrt(s) ~ 4.26 GeV
data_4360  = DatasetManager.real_data.find("703_4360")      # sqrt(s) ~ 4.36 GeV
data_4420  = DatasetManager.real_data.find("703_4420")      # sqrt(s) ~ 4.42 GeV

incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")

data_points  = [data_4230, data_4260, data_4360, data_4420]
incMC_points = [incMC_4230, incMC_4260, incMC_4360, incMC_4420]

### ------------------------- Decay cards (EvtGen) ------------------------- ###
# Mode 1 : eta_c -> K+ K- pi0
decay_card_KKpi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode 2 : eta_c -> K_S0 K+ pi-  (charge conjugate generated automatically for eta_c)
decay_card_KSKpi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0000 K_S0 K+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode 3 : eta_c -> 2(pi+ pi- pi0)
decay_card_2pipiPi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0000 pi+ pi- pi0 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode 4 : eta_c -> p+ anti-p-
decay_card_ppbar = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0000 p+ anti-p- PHSP;
  Enddecay

  End
DECAYCARD

# Inclusive mode : eta_c not reconstructed (only h_c -> gamma eta_c topology)
decay_card_inclusive = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 gamma eta_c PHSP;
  Enddecay

  End
DECAYCARD

### ------------------------- Exclusive MC (same signal over four energies) ------------------------- ###
exMC_KKpi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "hc_gamma_etac_KKpi0"
  config.events        = 100000
  config.decay_card    = decay_card_KKpi0
  config.cross_section = :default
end

exMC_KSKpi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "hc_gamma_etac_KSKpi"
  config.events        = 100000
  config.decay_card    = decay_card_KSKpi
  config.cross_section = :default
end

exMC_2pipiPi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "hc_gamma_etac_2pipiPi0"
  config.events        = 100000
  config.decay_card    = decay_card_2pipiPi0
  config.cross_section = :default
end

exMC_ppbar = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "hc_gamma_etac_ppbar"
  config.events        = 100000
  config.decay_card    = decay_card_ppbar
  config.cross_section = :default
end

### ============================ Event selection (BOSS) ============================ ###

# ---------------------------------------------------------------- Mode 1: K+ K- pi0
alg_name_1 = "HcToGammaEtaC_KKpi0"
alg_1 = Algorithm.new(alg_name_1)
alg_1.set_header(["#{alg_name_1}Alg/#{alg_name_1}.h"])
     .set_constant({"ECMS" => [:double, 4.230]})   # ECMS runs over 4.23/4.26/4.36/4.42 GeV
     .set_alias({"std::vector<double>" => "Vdouble"})
     .note(:hc_recoil_tag, "h_c tagged via recoil mass against pi+pi- in [3.515, 3.535] GeV/c^2; the pi+pi- tag pair is reconstructed in addition to the eta_c daughters")

sel_1 = Selection.new
  .select_track {
    cos_theta 0.93
    Vr        1.0
    Vz        10.0
    nChrp     "==2"      # pi+ (tag) + K+ : 2 positive tracks
    nChrn     "==2"      # pi- (tag) + K- : 2 negative tracks
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"    # 2 gamma from pi0 + 1 gamma from h_c
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]      # K+/-  vs pi
  }
  .remove([:kp <= :chrgp, :km <= :chrgn]) # strip identified kaons from charged lists
  .assign({:chrgp => :pip, :chrgn => :pim}) # remaining tracks are the tag pions
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # pi0 mass-constrained fit
    chi2_cut 25
    npi0  ">=1"
  }
  .kinematic_fit([:pip, :pim, :kp, :km, :pi0, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :kp, :km, :pi0).within(3.515, 3.535)  # h_c (recoil pi+pi-) window
    chi2_cut 25
  }

alg_1.with_decay_card(decay_card_KKpi0).apply(sel_1)
alg_1.execute_on(data_points + incMC_points + exMC_KKpi0)

# ---------------------------------------------------------------- Mode 2: K_S0 K+/- pi-/+ 
alg_name_2 = "HcToGammaEtaC_KSKpi"
alg_2 = Algorithm.new(alg_name_2)
alg_2.set_header(["#{alg_name_2}Alg/#{alg_name_2}.h"])
     .set_constant({"ECMS" => [:double, 4.230]})
     .set_alias({"std::vector<double>" => "Vdouble"})
     .note(:hc_recoil_tag, "h_c tagged via recoil mass against pi+pi- in [3.515, 3.535] GeV/c^2 (tag pair in addition to eta_c daughters)")

sel_2 = Selection.new
  .select_track {
    cos_theta 0.93
    Vr        1.0
    Vz        10.0
    nChrp     "==3"      # pi+ (tag) + K_S0 pi+ + K+ (eta_c)
    nChrn     "==3"      # pi- (tag) + K_S0 pi- + pi- (eta_c)
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"    # only the gamma from h_c -> gamma eta_c
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]      # K+/- vs pi
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:pip, :pim]) {                 # K_S0 -> pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:pip, :pim, :kp, :km, :K_S0, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 45
  }

alg_2.with_decay_card(decay_card_KSKpi).apply(sel_2)
alg_2.execute_on(data_points + incMC_points + exMC_KSKpi)

# ---------------------------------------------------------------- Mode 3: 2(pi+ pi- pi0)
alg_name_3 = "HcToGammaEtaC_2pipiPi0"
alg_3 = Algorithm.new(alg_name_3)
alg_3.set_header(["#{alg_name_3}Alg/#{alg_name_3}.h"])
     .set_constant({"ECMS" => [:double, 4.230]})
     .set_alias({"std::vector<double>" => "Vdouble"})
     .note(:hc_recoil_tag, "h_c tagged via recoil mass against pi+pi- in [3.515, 3.535] GeV/c^2; recoil window applied in ROOT where tag/eta_c pions overlap")

sel_3 = Selection.new
  .select_track {
    cos_theta 0.93
    Vr        1.0
    Vz        10.0
    nChrp     "==3"      # pi+ (tag) + 2 pi+ (eta_c)
    nChrn     "==3"      # pi- (tag) + 2 pi- (eta_c)
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=5"    # 4 gamma from 2 pi0 + 1 gamma from h_c
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon]      # pi+/- vs K
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # pi0 mass-constrained fit
    chi2_cut 25
    npi0  ">=1"
  }
  .kinematic_fit([:pip, :pim, :pip, :pip, :pim, :pim, :pi0, :pi0, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 35
  }

alg_3.with_decay_card(decay_card_2pipiPi0).apply(sel_3)
alg_3.execute_on(data_points + incMC_points + exMC_2pipiPi0)

# ---------------------------------------------------------------- Mode 4: p+ pbar-
alg_name_4 = "HcToGammaEtaC_ppbar"
alg_4 = Algorithm.new(alg_name_4)
alg_4.set_header(["#{alg_name_4}Alg/#{alg_name_4}.h"])
     .set_constant({"ECMS" => [:double, 4.230]})
     .set_alias({"std::vector<double>" => "Vdouble"})
     .note(:hc_recoil_tag, "h_c tagged via recoil mass against pi+pi- in [3.515, 3.535] GeV/c^2; tag pi+pi- reconstructed in addition to the p pbar pair")

sel_4 = Selection.new
  .select_track {
    cos_theta 0.93
    Vr        1.0
    Vz        10.0
    nChrp     "==2"      # pi+ (tag) + p
    nChrn     "==2"      # pi- (tag) + pbar
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"    # only the gamma from h_c -> gamma eta_c
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:pion]     # p / pbar vs pi
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim}) # remaining tracks are the tag pions
  .kinematic_fit([:pip, :pim, :prp, :prm, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :prp, :prm).within(3.515, 3.535)   # h_c (recoil pi+pi-) window
    chi2_cut 40
  }

alg_4.with_decay_card(decay_card_ppbar).apply(sel_4)
alg_4.execute_on(data_points + incMC_points + exMC_ppbar)

# ---------------------------------------------------------------- Inclusive (recoil) mode
alg_name_incl = "HcToGammaEtaC_Inclusive"
alg_incl = Algorithm.new(alg_name_incl)
alg_incl.set_header(["#{alg_name_incl}Alg/#{alg_name_incl}.h"])
        .set_constant({"ECMS" => [:double, 4.230]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:hc_recoil_tag, "h_c tagged via recoil mass against pi+pi- in [3.515, 3.535] GeV/c^2 (recoil_mass_of not expressible in the fit DSL; window applied in ROOT)")
        .note(:etac_recoil_tag, "eta_c tagged via recoil mass against pi+pi-gamma in [2.52, 3.4] GeV/c^2")

sel_incl = Selection.new
  .select_track {
    cos_theta 0.93
    Vr        1.0
    Vz        10.0
    nChrp     ">=1"      # at least 2 charged tracks (assumed pions)
    nChrn     ">=1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon]      # charged tracks assumed to be pions
  }
  # Partial reconstruction: reconstruct pi+ pi- gamma, infer the eta_c as recoil
  # recIDs : 1 = pi+, 2 = pi-, 4 = gamma(h_c)
  .partial_rec([1, 2, 4]) {
    require_recoil_mass 2.52, 3.4         # recoil against pi+pi-gamma  -> eta_c window
  }

alg_incl.with_decay_card(decay_card_inclusive).apply(sel_incl)
alg_incl.execute_on(data_points + incMC_points)