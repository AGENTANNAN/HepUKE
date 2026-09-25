### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # psi(3686) inclusive MC
cont_data  = DatasetManager.real_data.find("709_3650")       # 3.65 GeV continuum data (non-resonant bkg)

# --- Decay cards (EvtGen format) ---
# Signal: psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K^+- pi^-+ pi+ pi-, K_S0 -> pi+ pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c(2S) PHSP;
    Enddecay
    Decay eta_c(2S)
    1.000 K_S0 K+ pi- pi+ pi- PHSP;
    Enddecay
    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Background: psi(3686) -> gamma K_S0 K^+- pi^-+ pi+ pi- (phase space), K_S0 -> pi+ pi-
decay_card_ps = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma K_S0 K+ pi- pi+ pi- PHSP;
    Enddecay
    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# Background: psi(3686) -> pi0 K_S0 K^+- pi^-+ pi+ pi-, K_S0 -> pi+ pi-
decay_card_pi0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 K_S0 K+ pi- pi+ pi- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

# --- Exclusive MC samples (1M events each) ---
exmc_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gam_etac2s_ksk3pi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exmc_ps = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gam_ksk3pi_phsp"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_ps
  config.cross_section   = :default
end

exmc_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0_ksk3pi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_pi0
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name  = "PsiPToGamEtaC2S"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
  .select_track {                       # charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==3"                     # exactly three positive tracks
    nChrn     "==3"                     # exactly three negative tracks
    nNet      "==0"                     # net charge zero
  }
  .select_photon {                      # photon selection
    tdc_emc_start     0                 # EMC timing window 0-14
    tdc_emc_end       14
    energyThreshold_b 0.025            # E > 25 MeV in barrel
    energyThreshold_e 0.025            # E > 25 MeV in endcap
    nGam              ">=1"            # at least one photon
  }
  # Treat all charged tracks as pions to search all opposite-charge pion pairs
  .assign({:chrgp => :pip, :chrgn => :pim})
  # K_S0 -> pi+ pi- via secondary vertex fit, keeping best mass-difference combination
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list   # daughters no longer reusable
  }
  # Remaining four tracks: combinatorial PID (dE/dx + TOF), K/pi hypotheses only
  .pid(method: :chi2_sum) {
    chi_min_cut 4                      # chi2_min cut 4
    identify :kaon, :pion              # only K and pi hypotheses
  }
  # Nominal 4C kinematic fit: gamma K_S0 K+ pi- pi+ pi-
  .kinematic_fit([:gamma, :K_S0, :kp, :pim, :pip, :pim]) {
    nominal                            # flag as the nominal fit
    constrain_four_momentum            # 4C energy-momentum constraint
    chi2_cut 200                       # loose chi2 < 200 (tight cut done in ROOT)
  }

# BOSS-side requirement not expressible by a dedicated DSL method
algorithm.note(:ks_mass_window, "K_S0 candidates required to satisfy |M(pi+pi-) - M_K_S0| < 10 MeV/c2; the secondary vertex fit keeps the pair with the smallest mass difference, and the +-10 MeV/c2 window is the candidate-acceptance guard on the resulting K_S0")

algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = algorithm.execute_on([psip_data, psip_incMC, cont_data, exmc_signal, exmc_ps, exmc_pi0])