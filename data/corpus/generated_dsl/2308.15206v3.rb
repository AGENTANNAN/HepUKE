# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data (BOSS 709, 3.686 GeV)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # Corresponding inclusive MC

# Decay card for the signal process ψ(3686) → K⁻ Λ Ξ̄⁺ (EvtGen syntax):
#   Ξ̄⁺ → Λ̄ π⁺ , Λ̄ → p̄ π⁺ , Λ → p π⁻
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000  K-  Lambda0  anti-Xi+        PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000  anti-Lambda0  pi+            PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+                 HypWK;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                      HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the full decay chain (1M events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_KLambdaXibar_excl"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipKLambdaXibar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # √s = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                 # Charged-track selection
    cos_theta  0.93               # |cosθ| < 0.93
    Vz         10.0               # |Vz| < 10 cm
    Vr         1.0                # Vr < 1 cm
    nChrp      ">=2"              # at least 2 positive tracks
    nChrn      ">=2"              # at least 2 negative tracks
    nTot       ">=4"              # at least 4 charged tracks
  }
  .pid(method: :probability) {    # PID by the probability method
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # exactly one (anti-)proton (from Λ̄)
    identify :kaon,   against: [:pion]          # exactly one K⁻ (from ψ decay)
    identify :pion,   against: [:kaon]          # exactly two π⁺
    nprm "==1"
    nkm  "==1"
    npip "==2"
  }
  .secondary_vertex_fit([:prm, :pip]) {   # Λ̄ → p̄ π⁺ secondary vertex
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:Lambda_bar, :pip]) {   # Ξ̄⁺ → Λ̄ π⁺ secondary vertex
    build_virtual_particle(:Xi_bar_plus).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # 4C kinematic fit to the K⁻ Ξ̄⁺ system, treating the Λ as missing
  .kinematic_fit([:km, :Xi_bar_plus]) {
    nominal
    miss_track_of :Lambda               # Λ not reconstructed — inferred from the recoil
    constrain_four_momentum             # 4C energy–momentum constraint
    invariant_mass_of(:Xi_bar_plus).within(1.315, 1.330)   # Ξ̄⁺ mass window [1.315,1.330] GeV/c²
    chi2_cut 200                        # χ² < 200 (loose; tight cut applied in ROOT)
  }

# BOSS-side selections that have no dedicated DSL construct
my_algorithm
  .note(:mass_window, "Λ̄(p̄π⁺) invariant-mass window [1.110, 1.121] GeV/c² applied at the Λ̄ secondary-vertex reconstruction level")
  .note(:decay_length, "Ξ̄⁺ decay length > 0.5 cm required at the Ξ̄⁺ secondary-vertex reconstruction level")
  .note(:track_cuts, "Ξ̄⁺ daughter tracks use looser |Vz| < 15 cm and Vr < 10 cm cuts instead of the prompt-track |Vz| < 10 cm, Vr < 1 cm")
  .note(:recoil_mass, "prompt-Λ selection: recoil mass RM(K⁻Ξ̄⁺) required in [1.080, 1.140] GeV/c²; the Λ is reconstructed from the recoil against the K⁻Ξ̄⁺ system")
  .note(:charge_conjugate, "charge-conjugate channel ψ(3686) → K⁺ Λ̄ Ξ⁻ is analysed symmetrically")

# Attach the decay card and render the full event-selection chain
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Run on real data, inclusive MC, and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])