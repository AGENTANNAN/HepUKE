### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # ψ(2S) real data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # ψ(2S) inclusive MC
cont_data  = DatasetManager.real_data.find("709_3650")       # 3.65 GeV continuum real data
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")    # 3.65 GeV continuum inclusive MC

# --- Decay cards (EvtGen) ---
# Mode I: ψ(2S) → γ η_c(2S), η_c(2S) → K_S0 K+ π− (K_S0 → π+π−)
decay_card_kskpi = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c(2S) PHSP;
    Enddecay

    Decay eta_c(2S)
    1.000 K_S0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: ψ(2S) → γ η_c(2S), η_c(2S) → K+ K− π0 (π0 → γγ)
decay_card_kkpi0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c(2S) PHSP;
    Enddecay

    Decay eta_c(2S)
    1.000 K+ K- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC: 1,000,000 events per decay mode ---
exMC_kskpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_etac2S_KSKPi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_kskpi
  config.cross_section   = :default
end

exMC_kkpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_etac2S_KKPi0"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_kkpi0
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ================= Mode I: η_c(2S) → K_S0 K± π∓ =================
alg_name_kskpi = "Etac2SToKSKPi"
alg_kskpi = Algorithm.new(alg_name_kskpi)
alg_kskpi.set_header(["#{alg_name_kskpi}Alg/#{alg_name_kskpi}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .note(:ks0_selection_cuts,
               "K_S0 secondary-vertex candidates required to satisfy |M(pi+pi-) - m_K_S0| < 7 MeV and
                decay length > 2x the vertex resolution; these cuts act on the rebuilt K_S0 candidate
                and have no dedicated DSL method")
         .note(:modified_kinematic_fit,
               "modified 3C kinematic fit to gamma K_S0 K+ pi- : the radiative-photon energy is floated
                (left unconstrained) and the K_S0 K+ pi- system is constrained to the four-momentum of
                the recoil against the photon; loose chi2_m3C < 200 applied in BOSS, tight chi2_m3C < 20
                applied in the ROOT analysis")
         .note(:continuum_subtraction,
               "3.65 GeV continuum data (and its inclusive MC) used to subtract the continuum background")
         .note(:fsr_correction,
               "final-state radiation (FSR) correction applied to the signal")

sel_kskpi = Selection.new
  .select_track {
    cos_theta 0.93    # |cos(theta)| < 0.93
    Vz        10.0    # |Vz| < 10 cm
    Vr        1.0     # Vr < 1 cm
    nChrp     "==2"   # two positive charged tracks
    nChrn     "==2"   # two negative charged tracks
    nNet      "==0"   # net charge zero
  }
  .select_photon {
    energyThreshold_b 0.025   # barrel energy threshold 25 MeV
    energyThreshold_e 0.050   # endcap energy threshold 50 MeV
    tdc_emc_start     0       # TDC window start
    tdc_emc_end       14      # TDC window end
    angle_to_track    10.0    # angle to nearest charged track > 10 deg
    nGam              ">=1"   # at least one photon
  }
  .pid(method: :probability) {
    prob_cut 0.001                              # PID probability > 0.001
    identify :kaon, against: [:pion, :proton]   # separate K from pi/p (K+ and K-)
    nkp ">=1"                                   # at least one K+
    nkm ">=1"                                   # at least one K-
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])       # remove identified kaons from the charged-track lists
  .assign({:chrgp => :pip, :chrgn => :pim})     # treat remaining charged tracks as pi+ / pi-
  .secondary_vertex_fit([:pip, :pim]) {         # fit pi+pi- to a common vertex → K_S0
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :K_S0, :kp, :pim]) {  # modified 3C fit to gamma K_S0 K+ pi-
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_kskpi.with_decay_card(decay_card_kskpi).apply(sel_kskpi)


# ================= Mode II: η_c(2S) → K+ K− π0 =================
alg_name_kkpi0 = "Etac2SToKKPi0"
alg_kkpi0 = Algorithm.new(alg_name_kkpi0)
alg_kkpi0.set_header(["#{alg_name_kkpi0}Alg/#{alg_name_kkpi0}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .note(:modified_kinematic_fit,
               "modified 4C kinematic fit to gamma K+ K- pi0 : the radiative-photon energy is floated
                (left unconstrained) and the K+ K- pi0 system is constrained to the four-momentum of
                the recoil against the photon; loose chi2_m4C < 200 applied in BOSS, tight chi2_m4C < 15
                applied in the ROOT analysis")
         .note(:continuum_subtraction,
               "3.65 GeV continuum data (and its inclusive MC) used to subtract the continuum background")
         .note(:fsr_correction,
               "final-state radiation (FSR) correction applied to the signal")

sel_kkpi0 = Selection.new
  .select_track {
    cos_theta 0.93    # |cos(theta)| < 0.93
    Vz        10.0    # |Vz| < 10 cm
    Vr        1.0     # Vr < 1 cm
    nChrp     "==1"   # one positive charged track
    nChrn     "==1"   # one negative charged track
    nNet      "==0"   # net charge zero
  }
  .select_photon {
    energyThreshold_b 0.025   # barrel energy threshold 25 MeV
    energyThreshold_e 0.050   # endcap energy threshold 50 MeV
    tdc_emc_start     0       # TDC window start
    tdc_emc_end       14      # TDC window end
    angle_to_track    10.0    # angle to nearest charged track > 10 deg
    nGam              ">=3"   # at least three photons (1 radiative + 2 from pi0)
  }
  .pid(method: :probability) {
    prob_cut 0.001                              # PID probability > 0.001
    identify :kaon, against: [:pion, :proton]   # separate K from pi/p (K+ and K-)
    nkp "==1"                                   # exactly one K+
    nkm "==1"                                   # exactly one K-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # rebuild pi0 from two photons (1C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 15
    npi0 ">=1"                                  # at least one pi0 candidate
  }
  .kinematic_fit([:gamma, :kp, :km, :pi0]) {    # modified 4C fit to gamma K+ K- pi0
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_kkpi0.with_decay_card(decay_card_kkpi0).apply(sel_kkpi0)


### Execute on datasets ###
root_files_kskpi = alg_kskpi.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_kskpi])
root_files_kkpi0 = alg_kkpi0.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_kkpi0])