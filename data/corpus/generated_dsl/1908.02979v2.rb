# ============================================================
# Dataset preparation
# ============================================================
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(2S) 3.686 GeV real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # matching inclusive MC

# Decay card for the RADIATIVE mode
#   psi(2S) -> gamma chi_cJ (E1 transitions, J = 0,1,2)
#   chi_cJ  -> anti-p- K*+ Lambda0
#   K*+     -> K+ pi0 ;  pi0 -> gamma gamma ;  Lambda0 -> p+ pi-
decay_card_radiative = <<~DECAYCARD
    Decay psi(2S)
    0.3333 gamma chi_c0 P2GC0;
    0.3333 gamma chi_c1 P2GC1;
    0.3334 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c0
    1.000 anti-p- K*+ Lambda0 PHSP;
    Enddecay

    Decay chi_c1
    1.000 anti-p- K*+ Lambda0 PHSP;
    Enddecay

    Decay chi_c2
    1.000 anti-p- K*+ Lambda0 PHSP;
    Enddecay

    Decay K*+
    1.000 K+ pi0 VSS;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    End
DECAYCARD

# Decay card for the DIRECT mode
#   psi(2S) -> anti-p- K*+ Lambda0
#   K*+     -> K+ pi0 ;  pi0 -> gamma gamma ;  Lambda0 -> p+ pi-
decay_card_direct = <<~DECAYCARD
    Decay psi(2S)
    1.000 anti-p- K*+ Lambda0 PHSP;
    Enddecay

    Decay K*+
    1.000 K+ pi0 VSS;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for each of the two signal modes
exMC_radiative = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_to_gamma_chicJ_to_antip_Kstar_Lambda"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_radiative
  config.cross_section   = :default
end

exMC_direct = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_to_antip_Kstar_Lambda"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_direct
  config.cross_section   = :default
end

# ============================================================
# Event selection (BOSS)
# ============================================================

# ---- Radiative channel: psi(2S) -> gamma chi_cJ -> gamma anti-p K+ Lambda pi0 ----
alg_name_rad = "PsiPToGammaChiCJ"
alg_rad = Algorithm.new(alg_name_rad)
alg_rad.set_header(["#{alg_name_rad}Alg/#{alg_name_rad}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})          # 3.686 GeV center-of-mass energy
       .set_alias({"std::vector<double>" => "Vdouble"})
alg_rad.note(:pid_correction_method,
  "The low-momentum pion from Lambda -> p+ pi- falls below the TOF reachable momentum, so only
   dE/dx information is used for its identification; TOF is not applied to this daughter pion.")

sel_rad = Selection.new
  .select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93
    Vz        100.0     # |Vz| < 100 cm
    Vr        10.0      # Vr < 10 cm
    nChrp     ">=2"     # at least two positive tracks
    nChrn     ">=2"     # at least two negative tracks
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025   # barrel energy threshold 25 MeV
    energyThreshold_e 0.050   # endcap energy threshold 50 MeV
    angle_to_track    5.0     # > 5 deg from any charged track
    nGam              ">=3"   # at least three photons (E1 photon + two from pi0)
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]    # p+ and anti-p, against K and pi
    identify :kaon,   against: [:pion, :proton]  # K+, against pi and p
    nprp ">=1"     # at least one proton
    nprm ">=1"     # at least one anti-proton
    nkp  ">=1"     # at least one K+
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn, :kp <= :chrgp])  # drop used tracks from charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})                 # remaining tracks as pi+/pi-
  .select_isolated_photon {                 # isolated-photon requirement (10 deg from primary (anti)protons)
    angle_to_prp_track 10.0
    angle_to_prm_track 10.0
    nGam ">=3"
  }
  .secondary_vertex_fit([:prp, :pim]) {     # Lambda -> p+ pi- secondary vertex fit
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {                  # 1-C pi0 reconstruction: gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"     # at least one pi0 candidate
  }
  .kinematic_fit([:gamma, :prm, :kp, :Lambda, :pi0]) {       # nominal 5C fit (radiative hypothesis)
    nominal
    constrain_four_momentum
    chi2_cut 70
  }
  .kinematic_fit([:prm, :kp, :Lambda, :pi0]) {               # competing 4C (direct) hypothesis
    constrain_four_momentum                                  # store chi2 only -> ROOT-level veto
  }

alg_rad.with_decay_card(decay_card_radiative).apply(sel_rad)
root_files_rad = alg_rad.execute_on([psip_data, psip_incMC, exMC_radiative])

# ---- Direct channel: psi(2S) -> anti-p K+ Lambda pi0 ----
alg_name_dir = "PsiPToAntiPKStarLambda"
alg_dir = Algorithm.new(alg_name_dir)
alg_dir.set_header(["#{alg_name_dir}Alg/#{alg_name_dir}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})
alg_dir.note(:pid_correction_method,
  "The low-momentum pion from Lambda -> p+ pi- falls below the TOF reachable momentum, so only
   dE/dx information is used for its identification; TOF is not applied to this daughter pion.")

sel_dir = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     ">=2"
    nChrn     ">=2"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    5.0
    nGam              ">=2"   # at least two photons (from pi0)
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion, :proton]
    nprp ">=1"
    nprm ">=1"
    nkp  ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn, :kp <= :chrgp])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .select_isolated_photon {
    angle_to_prp_track 10.0
    angle_to_prm_track 10.0
    nGam ">=2"
  }
  .secondary_vertex_fit([:prp, :pim]) {     # Lambda -> p+ pi- secondary vertex fit
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {                  # 1-C pi0 reconstruction: gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:prm, :kp, :Lambda, :pi0]) {               # nominal 4C fit (direct hypothesis)
    nominal
    constrain_four_momentum
    chi2_cut 100
  }
  .kinematic_fit([:gamma, :prm, :kp, :Lambda, :pi0]) {       # competing 5C (radiative) hypothesis
    constrain_four_momentum                                  # store chi2 only -> ROOT-level veto
  }

alg_dir.with_decay_card(decay_card_direct).apply(sel_dir)
root_files_dir = alg_dir.execute_on([psip_data, psip_incMC, exMC_direct])