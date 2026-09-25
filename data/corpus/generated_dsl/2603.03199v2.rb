# ============================================================================
# BOSS part (dataset preparation + event selection) for the search of a
# massless BSM particle X produced in J/psi(3097) decays:
#
#   J/psi -> anti-Xi0 ( -> anti-Lambda0 pi0 , anti-Lambda0 -> anti-p- pi+ ,
#                         pi0 -> gamma gamma )            <-- TAG side
#            Xi0      ( -> Lambda0 + invisible(X) , Lambda0 -> p+ pi- )
#
# Final state: p pbar pi+ pi- pi0 + invisible
# ============================================================================

### Datasets ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi(3097) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

### Decay cards ###
# Full signal chain: Xi0 -> Lambda + X (X invisible/massless), tagged by anti-Xi0
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000  Xi0  anti-Xi0  PHSP;
  Enddecay

  Decay Xi0
  1.0000  Lambda0  X  PHSP;
  Enddecay

  Decay anti-Xi0
  1.0000  anti-Lambda0  pi0  PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-  HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000  anti-p-  pi+  HypWK;
  Enddecay

  Decay pi0
  1.0000  gamma  gamma  PHSP;
  Enddecay

  End
DECAYCARD

# Single-tag sample: anti-Xi0 -> anti-Lambda0 pi0 (other side left unspecified)
decay_card_singletag = <<~DECAYCARD
  Decay J/psi
  1.0000  anti-Xi0  X  PHSP;
  Enddecay

  Decay anti-Xi0
  1.0000  anti-Lambda0  pi0  PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000  anti-p-  pi+  HypWK;
  Enddecay

  Decay pi0
  1.0000  gamma  gamma  PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC samples (500k events each) ###
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_jpsi_xi0_invisible_antixi0_tag"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exMC_singletag = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exMC_jpsi_antixi0_singletag"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_singletag
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiXiInvisible"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    .select_track {                       # charged track selection
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz        20.0                    # |Vz| < 20 cm
        nChrp     ">=2"                   # >= 2 positive tracks (p and pi+)
        nChrn     ">=2"                   # >= 2 negative tracks (pbar and pi-)
        nNet      "==0"                   # net charge zero
    }
    .select_photon {                      # photon selection
        tdc_emc_start     0               # TDC window [0, 14]
        tdc_emc_end       14
        angle_to_track    10.0            # >= 10 deg from nearest (prompt) charged track
        energyThreshold_b 0.025           # 25 MeV in the barrel
        energyThreshold_e 0.050           # 50 MeV in the endcap
        nGam              ">=2"           # >= 2 photons
    }
    .pid(method: :probability) {          # PID: probability method, p / pi / K separation
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # p+ and pbar
        identify :pion,   against: [:kaon, :proton] # pi+ and pi-
        nprp ">=1"                        # >= 1 proton
        nprm ">=1"                        # >= 1 anti-proton
        npip ">=1"                        # >= 1 pi+
        npim ">=1"                        # >= 1 pi-
    }
    .select_isolated_photon {             # require >= 20 deg from anti-proton (pbar) tracks
        angle_to_prm_track 20.0
        nGam ">=2"
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C mass-constrained pi0 -> gamma gamma
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
    }
    .secondary_vertex_fit([:prm, :pip]) {       # anti-Lambda0 -> pbar pi+  (minimum mass difference)
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prp, :pim]) {       # Lambda0 -> p pi-       (minimum mass difference)
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # ----------------------------------------------------------------------
    # Nominal 4-momentum-conserving fit: J/psi -> p pbar pi+ pi- pi0 + invisible
    # (the invisible BSM particle X is treated as a missing track)
    # ----------------------------------------------------------------------
    .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :X]) {
        nominal
        constrain_four_momentum
        miss_track_of :X
        chi2_cut 200                      # loose cut only; tight value optimised in ROOT
    }
    # Competing 2C hypothesis: the invisible recoil is instead a pi0
    # (SM Xi0 -> Lambda pi0). Non-nominal -> chi2 stored for the ROOT-level veto
    # chi2(invisible-pi0) > chi2(massless).
    .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :pi0]) {
        constrain_four_momentum
    }
    # Extra-photon hypotheses (5C / 6C), non-nominal -> chi2 stored; the ROOT-level
    # veto requires min chi2 > 1000.
    .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :gamma, :X]) {
        constrain_four_momentum
        miss_track_of :X
    }
    .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :gamma, :gamma, :X]) {
        constrain_four_momentum
        miss_track_of :X
    }

# BOSS-side procedures that cannot be expressed in formal DSL syntax
my_algorithm
    .note(:background_veto,
          "competing 2C fits: the invisible-recoil hypothesis is compared with a pi0 "
          "hypothesis and with an extra-photon (5C/6C) hypothesis. The requirement "
          "chi2(2C) < 8.5 and chi2(invisible-pi0) > chi2(massless), as well as the "
          "extra-photon veto min chi2 > 1000, are applied at the ROOT level on the "
          "stored chi2 values.")
    .note(:efficiency_curve,
          "data-driven E_extra shape correction applied to the extra-photon energy "
          "spectrum (control sample derived) before the 5C/6C extra-photon veto")
    .with_decay_card(decay_card_signal)
    .apply(event_selection)

# Run on real data, inclusive MC and both exclusive MC samples
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_singletag])