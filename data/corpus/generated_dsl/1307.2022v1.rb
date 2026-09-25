# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
# psi(3770) energy point (sqrt(s) = 3.773 GeV)
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
# Continuum energy point (sqrt(s) = 3.650 GeV)
data_3650  = DatasetManager.real_data.find("709_3650")
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")

# Decay card for the Bhabha signal e+e- -> (gamma) e+ e-.
# The description asks for the Babayaga generator; Babayaga cannot be expressed in an
# EvtGen decay card, so the documented KKMC + psi(4260) fallback top-mother is used
# (the generator settings are preserved in the note block below).
decay_card_bhabha = <<~DECAYCARD
    Decay psi(4260)
    1.0 e+ e- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

# 400k-event exclusive Bhabha MC produced at BOTH energy points (3.650 and 3.773 GeV).
# create_exclusive_mc_for returns one ExclusiveMC per related dataset (same card and
# event count, differing only by the associated real dataset).
exMC_bhabha = DatasetManager.create_exclusive_mc_for([data_3650, data_3773]) do |config|
  config.sample_name   = "exmc_bhabha"
  config.events        = 400_000
  config.decay_card    = decay_card_bhabha
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name   = "Bhabha"
bhabha_alg = Algorithm.new(alg_name)
bhabha_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
          .set_constant({"ECMS" => [:double, 3.7730]})
          .set_alias({"std::vector<double>" => "Vdouble"})

# Both energy points share this single selection chain.
event_selection = Selection.new
  .select_track {                 # exactly one e+ and one e- from a common vertex
      cos_theta 0.80              # |cos(theta)| < 0.80
      Vz        5.0               # |Vz| < 5 cm
      Vr        10.0              # |Vr| < 10 mm (1 cm) at closest approach
      nChrp     "==1"             # exactly one positively charged track
      nChrn     "==1"             # exactly one negatively charged track
      nNet      "==0"             # net charge zero
  }
  .pid(method: :probability) {    # lepton identification with the probability method
      prob_cut 0.001              # PID probability > 0.001
      # track with p > 1.0 GeV -> lepton; lepton with EMC energy > 0.6 GeV -> electron,
      # otherwise muon
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      nlp "==1"                   # exactly one positive lepton (e+)
      nlm "==1"                   # exactly one negative lepton (e-)
  }
  .kinematic_fit([:lp, :lm]) {    # nominal 4C fit of the e+e- system
      nominal                     # mandatory nominal fit of the selection chain
      constrain_four_momentum     # constrain the e+e- four-momentum to the CMS energy
      chi2_cut 200                # loose cut; the tight cut is chosen later in ROOT
  }

# BOSS-side procedures that have no formal DSL representation
bhabha_alg
  .note(:emc_energy_cut, "both the positron and the electron candidate are required
    to have EMC cluster energy > 1.0 GeV in order to reject mu+mu- events; the DSL
    exposes no accessor for the EMC energy of a selected lepton track, so the cut is
    applied in the generated BOSS code right after the lepton identification step")
  .note(:cosmic_ray_veto, "cosmic-ray events are rejected by requiring that the two
    track momenta |p1| and |p2| are not both larger than E_beam + 0.15 GeV with
    E_beam = sqrt(s)/2; this event-level two-track condition has no DSL form and is
    applied in the generated BOSS code")
  .note(:background_veto, "radiative (gamma)J/psi, (gamma)psi(3686) and
    psi(3770) -> (gamma) J/psi X backgrounds are suppressed by requiring
    |p1| + |p2| > 0.9 * sqrt(s); implemented as an event-level cut in the generated
    BOSS code")
  .note(:generator_settings, "the exclusive Bhabha MC is generated with the Babayaga
    generator (ISR/FSR included) with a generator-level cut |cos(theta)| < 0.83; the
    EvtGen decay card can only carry the KKMC psi(4260) fallback, so the Babayaga
    settings must be injected into the generator jobOptions")
  .note(:beam_energy_per_point, "one algorithm serves two energy points (3.650 GeV
    continuum and 3.773 GeV psi(3770)); ECMS is set to 3.7730 GeV, so the 3.650 GeV
    sample must be processed with ECMS (and hence P4_cms) overridden to 3.650 GeV")
  .note(:trigger_efficiency, "the trigger efficiency for e+e- -> (gamma)e+e- is taken
    to be 100%, i.e. no trigger correction is applied to the selected events")
  .with_decay_card(decay_card_bhabha)
  .apply(event_selection)

# Execute on both energy points (real data + inclusive MC) and on the exclusive
# Bhabha MC generated at both energies.
root_files = bhabha_alg.execute_on([data_3773, incMC_3773, data_3650, incMC_3650] + exMC_bhabha)