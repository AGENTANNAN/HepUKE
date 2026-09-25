### Dataset description ###
data_3097  = DatasetManager.real_data.find("708_3097")        # J/psi(3097) real data
incMC_3097 = DatasetManager.inclusive_mc.find("708_3097")      # Corresponding inclusive MC

### Decay card for the full signal chain (EvtGen format) ###
# J/psi -> Sigma+ anti-Sigma-,
#   Sigma+     -> Lambda e+ nu_e   (Lambda -> p+ pi-),
#   anti-Sigma- -> anti-p- pi0     (pi0 -> gamma gamma)
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 Sigma+ anti-Sigma- PHSP;
  Enddecay

  Decay Sigma+
  1.0000 Lambda e+ nu_e PHSP;
  Enddecay

  Decay Lambda
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Sigma-
  1.0000 anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC for the full decay chain (1,000,000 events) ###
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_SigmaSigma_LambdaeNu"
  config.related_dataset = data_3097
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "SigmaSigmaDTag"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV

selection = Selection.new
  # Exactly four charged tracks (p, e+ ; anti-p, pi-)
  .select_track {
    cos_theta 0.93       # |cos(theta)| < 0.93
    nChrp     "==2"      # exactly 2 positive tracks  (p, e+)
    nChrn     "==2"      # exactly 2 negative tracks  (anti-p, pi-)
    nNet      "==0"      # net charge zero
  }
  # At least two photons (from pi0 -> gamma gamma)
  .select_photon {
    tdc_emc_start     0      # TDC window [0, 14]
    tdc_emc_end       14
    angle_to_track    10.0   # angle to nearest charged track > 10 deg
    energyThreshold_b 0.025  # barrel threshold 25 MeV
    energyThreshold_e 0.050  # endcap threshold 50 MeV
    nGam              ">=2"
  }
  # pi0 reconstruction: 1C kinematic fit of gamma gamma to the pi0 nominal mass
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0     ">=1"
  }
  # PID (probability method)
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]        # L(p) > L(K), L(p) > L(pi), L(p) > 0.001
    identify :pion,   against: [:kaon, :electron]    # L(pi) > L(K), L(pi) > L(e)
    # high-momentum e+ (leptons); the EMC+TOF likelihood ratio is not DSL-tunable (see note)
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    nprp "==1"
    nprm "==1"
    npim "==1"
    nlp  "==1"
  }
  # Lambda -> p pi- secondary vertex (anti-p and e+ left untouched for the tag/signal sides)
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Final kinematic fit: tag (anti-p pi0) + signal (Lambda e+) with the neutrino missed
  .kinematic_fit([:prm, :pi0, :Lambda, :lp]) {
    nominal
    miss_track_of :nu_e          # neutrino treated as a missing particle
    constrain_four_momentum
    chi2_cut 200                 # loose BOSS-side cut; tight cut in ROOT
  }

# BOSS-side procedures that have no dedicated DSL primitive
alg.note(:tag_selection_windows,
         "ST/DT tag windows for the anti-Sigma- -> anti-p pi0 tag are applied as "
       + "selection: ST deltaE_tag in [-48, 52] MeV (keep the candidate with minimum "
       + "|deltaE_tag|) and DT |M_BC^tag - m_Sigma-| < 0.049 GeV/c^2.")
   .note(:missing_mass_construction,
         "The neutrino is treated as missing; M_miss^2 is built from "
       + "E_miss = E_beam - E_Lambda - E_e+ and p_miss = |p_Sigma+ - p_Lambda - p_e+|, "
       + "with the Sigma+ momentum inferred from the tag as "
       + "-p_hat_Sigma- * sqrt(E_beam^2 - m_Sigma+^2).")
   .note(:background_veto,
         "Tag-side background suppression: angle(anti-p, pi0) > 170 deg and "
       + "p(anti-p) in [0.16, 0.21] GeV/c in the anti-Sigma- rest frame; "
       + "M_recoil(Sigma, Lambda) + M(p pi-) - m_Lambda > -60 MeV/c^2.")
   .note(:misid_veto,
         "Veto M_Lambda(pi+ -> e+) < 1.32 GeV/c^2 against e+ mis-identified as pi+.")
   .note(:vertex_quality,
         "Lambda secondary-vertex candidates additionally require "
       + "|M(p pi-) - m_Lambda| < 10 MeV/c^2 and a decay length larger than twice the "
       + "vertex resolution from the IP.")
   .note(:pid_correction_method,
         "The high-momentum e+ is selected with the combined EMC+TOF likelihood ratio "
       + "L'(e)/(L'(e)+L'(pi)+L'(K)) > 0.8; identify_high_momentum_leptons is used here "
       + "as the closest available DSL primitive.")

# Generate the algorithm for this decay card and run it
alg.with_decay_card(decay_card_signal).apply(selection)
root_files = alg.execute_on([data_3097, incMC_3097, exMC_signal])