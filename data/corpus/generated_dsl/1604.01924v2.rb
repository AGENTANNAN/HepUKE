# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# psi(3770) real data and the matching inclusive MC sample
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: e+e- -> mu+ mu- gamma via initial-state radiation at sqrt(s)=3.773 GeV.
# PHOKHARA models the ISR photon emission over the full photon-angle range; no
# intermediate charmonium (e.g. J/psi) is imposed at generator level.
decay_card_isr = <<~DECAYCARD
    Decay psi(3770)
    1.0000 mu+ mu- gamma PHOKHARA;
    Enddecay
    End
DECAYCARD

# Exclusive signal MC: 100k events of e+e- -> mu+ mu- gamma (PHOKHARA)
exMC_isr = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_mumugamma_isr"
  config.related_dataset = data_3773          # match the real 3.773 GeV sample
  config.events          = 100000
  config.decay_card      = decay_card_isr
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiISRMuMu"
jpsi_isr = Algorithm.new(alg_name)
jpsi_isr.set_header(["#{alg_name}Alg/#{alg_name}.h"])
        .set_constant({"ECMS" => [:double, 3.773]})  # sqrt(s) = 3.773 GeV
        # BOSS-side procedures that have no formal DSL construct:
        .note(:background_veto, "hadronic background is further rejected by a TMVA ANN; both muon candidate tracks must satisfy y_ANN < 0.3")
        .note(:track_pair_selection, "in three-track events the muon pair closest to the interaction point is used")
        .note(:pt_cut, "each selected lepton track must have transverse momentum p_t > 300 MeV/c")
        .note(:efficiency_curve, "selection efficiency taken from the full-range true-MC J/psi sample with per-track corrections, giving (32.04 +/- 0.09)%")

event_selection = Selection.new
event_selection
  .select_track {
      cos_theta  0.93    # |cos(theta)| < 0.93
      Vz         10.0    # |Vz| < 10 cm
      Vr         1.0     # Vr < 1 cm in the transverse plane
      nChrp      ">=1"   # at least one positive track
      nChrn      ">=1"   # at least one negative track
      nNet       "==0"   # net charge zero (>= 2 charged tracks in total)
  }
  .pid(method: :probability) {
      prob_cut 0.001     # PID probability > 0.001
      # High-momentum tracks (p > 1.0 GeV/c) treated as leptons; lepton classified
      # as electron if EMC energy > 0.6 GeV, otherwise as muon (P(mu) > P(e)).
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      nmup "==1"         # exactly one mu+
      nmum "==1"         # exactly one mu-
  }
  # Nominal 1C kinematic fit: e+e- -> mu+ mu- gamma with a missing massless photon
  # (four-momentum balance). The ISR photon is untagged and left missing.
  .kinematic_fit([:mup, :mum, :gamma]) {
      nominal
      miss_track_of :gamma        # photon is not measured (massless)
      constrain_four_momentum
      chi2_cut 10                 # chi2 < 10
  }
  # Competing hypothesis without the ISR photon; chi2 stored for a later
  # (ROOT-level) consistency check -- no chi2_cut, no nominal (Rule T2).
  .kinematic_fit([:mup, :mum]) {
      constrain_four_momentum
  }

# Generate the BOSS algorithm for the ISR signal process in the decay card.
jpsi_isr.with_decay_card(decay_card_isr).apply(event_selection)

# Execute on real data, inclusive MC, and the exclusive ISR signal MC.
root_files = jpsi_isr.execute_on([data_3773, incMC_3773, exMC_isr])