# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC

# Decay card: ψ(3686) → e+ μ−  (the charge-conjugate e− μ+ mode is handled analogously)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 e+ mu- PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the CLFV signal
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_emu"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipToEMu"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

# Pre-selection (tracks + photons) → lepton PID → 4C kinematic fit.
# The 4C fit is the final fit; the back-to-back (Δθ, Δφ), Bhabha veto and the
# |Σp|/√s, ΣE/√s signal-region cuts are applied after the fit and therefore
# belong to the ROOT analysis stage.
event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta 0.93        # |cosθ| < 0.93
                  Vz        10.0        # |Vz| < 10 cm
                  Vr        1.0         # Vr < 1 cm
                  nChrp     "==1"       # exactly one positive track
                  nChrn     "==1"       # exactly one negative track
                  nNet      "==0"       # net charge zero (exactly two tracks)
                }
               .select_photon {         # Photon selection — no photons allowed
                  tdc_emc_start     0   # TDC window [0, 700] ns
                  tdc_emc_end       14
                  angle_to_track    10.0  # angle to nearest charged track > 10°
                  energyThreshold_b 0.025 # > 25 MeV (barrel)
                  energyThreshold_e 0.050 # > 50 MeV (endcap)
                  nGam              "==0"
                }
               .pid(method: :probability) {   # Probability PID
                  prob_cut 0.001
                  # p > 1.0 GeV tracks treated as leptons, classified as electrons by the
                  # E/p-like discriminant (approximated here by the EMC energy threshold)
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.94
                  nlp "==1"   # one positive lepton
                  nlm "==1"   # one negative lepton
                }
               # 4C kinematic fit of the e± μ∓ pair (nominal, final fit)
               .kinematic_fit([:lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

# The E/p > 0.94 electron criterion used by the high-momentum lepton ID is not
# directly expressible — identify_high_momentum_leptons only exposes an EMC-energy
# threshold; the exact discriminant is preserved here for the systematic-uncertainty study.
my_algorithm.note(:electron_id_criterion,
  "high-momentum tracks (p > 1.0 GeV) are classified as electrons when E/p > 0.94; " \
  "the DSL's identify_high_momentum_leptons exposes only an EMC-energy threshold, so the " \
  "E/p discriminant itself is applied/validated outside the DSL")

# Generate the algorithm for the decay card and run on data, inclusive MC and signal MC
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])