# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # Corresponding inclusive MC at 3.773 GeV

# Decay card for the signal process (EvtGen format):
#   psi(3770) -> D0 D0bar,
#   D0    -> K_S0 K+ K-  (the untagged signal side),
#   K_S0  -> pi+ pi-,
#   D0bar -> K+ pi-      (the accompanying, non-reconstructed side).
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 K+ K- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 200k exclusive signal MC events for psi(3770) -> D0 D0bar
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKsKPi"
  config.related_dataset = data_3773        # match simulation conditions to the real data
  config.events          = 200000           # 200k events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "D0ToKsKPi"                       # untagged D0 -> K_S0 K+ K-
my_algorithm = Algorithm.new(alg_name)
my_algorithm
  .set_header(["#{alg_name}Alg/#{alg_name}.h"])
  .set_constant({"ECMS" => [:double, 3.773]}) # CMS energy 3.773 GeV
  # --- BOSS-side procedures that cannot be expressed in the DSL are captured below ---
  .note(:ks0_selection, "K_S0 candidates are required to have invariant mass in [0.487, 0.511] GeV/c^2, "
                        "secondary-vertex chi2 < 100 and flight-distance significance > 2 (untagged sample). "
                        "These criteria are not expressible in the secondary_vertex_fit block and are applied "
                        "at the ROOT level.")
  .note(:candidate_selection, "When several D0 -> K_S0 K+ K- candidates survive, the best one is chosen by the "
                              "smallest |E_D0 - E_beam|. The DSL kinematic_fit ranks combinations by chi2, so the "
                              "|E_D0 - E_beam| ranking is applied at the ROOT level.")
  .note(:track_selection_uniform, "The same track selection (|cos_theta| < 0.93, |Vz| < 10 cm, Vr < 1 cm) and PID "
                                  "are applied to all charged tracks, including the K_S0 daughters; the relaxed "
                                  "K_S0-specific tracking (Vz < 20 cm, no Vr, no PID for K_S0 pions) is intentionally "
                                  "not separately encoded.")

# Single untagged selection chain (no tagged / Dalitz selection).
event_selection = Selection.new
event_selection
  .select_track {                              # Charged track selection
    cos_theta 0.93                             # |cos(theta)| < 0.93
    Vz        10.0                             # |Vz| < 10 cm
    Vr        1.0                              # Vr < 1 cm
    nChrp     ">=2"                            # at least two positive tracks
    nChrn     ">=2"                            # at least two negative tracks
  }
  .pid(method: :probability) {                 # PID with the probability method (kaons favoured over pions)
    prob_cut 0.001                             # PID probability > 0.001
    identify :kaon, against: [:pion, :proton]  # K+ and K- (charge-conjugation shorthand)
    identify :pion, against: [:kaon, :proton]  # pi+ and pi-
    nkp  ">=1"                                 # at least one K+
    nkm  ">=1"                                 # at least one K-
    npip ">=1"                                 # at least one pi+
    npim ">=1"                                 # at least one pi-
  }
  .secondary_vertex_fit([:pip, :pim]) {        # Build K_S0 from pi+ pi- via a secondary-vertex fit
    build_virtual_particle(:K_S0).by_minimizing_mass_difference  # choose combination closest to the K_S0 mass
    remove_used_particle_from_candidate_list   # remove used pions from their lists
  }
  # 1C kinematic fit to K_S0 K+ K- constraining the invariant mass to the nominal D0 mass
  .kinematic_fit([:K_S0, :kp, :km]) {
    nominal
    invariant_mass_of(:K_S0, :kp, :km).constrain_to_nominal_mass_of(:D0)
    chi2_cut 20                                # chi2 < 20
  }

# Generate the complete algorithm for the process in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Run the single untagged selection chain on real data, inclusive MC and exclusive signal MC
root_files = my_algorithm.execute_on([data_3773, incMC_3773, exMC_signal])