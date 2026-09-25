# ============================================================
# J/psi(3097) -> gamma a, a -> gamma gamma   (purely neutral 3-photon final state)
# ============================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi(3097) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# Signal decay card: J/psi -> gamma a, a -> gamma gamma (intermediate pseudoscalar "a")
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma a PHSP;
    Enddecay

    Decay a
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process (1,000,000 events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_a_gg"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name  = "JpsiGammaAGG"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
    .select_track {                 # purely neutral final state: zero charged tracks
        nChrp "==0"                 # no positively charged tracks
        nChrn "==0"                 # no negatively charged tracks
        nNet  "==0"                 # net charge zero
    }
    .select_photon {                # good photons, barrel only, at least three
        tdc_emc_start     0         # EMC TDC window 0-14
        tdc_emc_end       14
        angle_to_track    10.0      # photon-track separation angle > 10 degrees
        energyThreshold_b 0.025     # barrel energy threshold E > 25 MeV (|cos(theta)| < 0.80)
        energyThreshold_e 1.0e9     # endcap reconstruction disabled
        nGam              ">=3"     # at least three photons
    }                               # no PID applied (no charged tracks in the final state)
    .kinematic_fit([:gamma, :gamma, :gamma]) {   # nominal 4C kinematic fit to the three photons
        nominal
        constrain_four_momentum
        chi2_cut 30                 # chi2(4C) < 30
        invariant_mass_of(:gamma, :gamma).out_of(0.11, 0.16)   # veto pi0 window
        invariant_mass_of(:gamma, :gamma).out_of(0.52, 0.56)   # veto eta window
        invariant_mass_of(:gamma, :gamma).out_of(0.92, 0.99)   # veto eta' window
        invariant_mass_of(:gamma, :gamma).out_of(2.92, 3.04)   # veto eta_c window
    }
    # competing hypotheses: 2-, 4- and 5-photon 4C fits. No chi2_cut and no nominal,
    # so their chi2 values are stored for the ROOT-level veto
    #   chi2(3gamma) < chi2(2gamma), chi2(4gamma), chi2(5gamma)
    .kinematic_fit([:gamma, :gamma]) {
        constrain_four_momentum
    }
    .kinematic_fit([:gamma, :gamma, :gamma, :gamma]) {
        constrain_four_momentum
    }
    .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma]) {
        constrain_four_momentum
    }

# BOSS-side background rejections that have no DSL construct
algorithm
    .note(:background_veto,
          "Delta E13 < 1.46 GeV and Delta E23 < 1.41 GeV for photons ordered by decreasing energy, and |Delta phi31| > 1 rad, applied after the nominal 4C fit; these ordering-dependent photon energy / azimuthal-angle rejections are not expressible in the DSL and are applied at ROOT level")
    .note(:photon_pair_timing,
          "photon-pair timing requirement |Delta T| < 500 ns on the EMC TDC; not expressible in the DSL")
    .with_decay_card(decay_card_signal)
    .apply(event_selection)

# Execute on real data, inclusive MC and signal exclusive MC
root_files = algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])