# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")      # ψ(3770) 3.773 GeV real data
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC sample

# Decay card for the signal process (EvtGen format, EvtGen particle names)
# ψ(3770) → D0 anti-D0 ; signal side D0 → K+K−π+π− ; tagged side anti-D0 → K+π− (charge conjugation implied)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K+ K- pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process: 1M events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKKPiPi_D0barToKPi"
  config.related_dataset = psi3770_data      # associated real dataset (for better simulation)
  config.events          = 1_000_000         # 1M events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) — tag-based analysis (TagAnalysis, not Selection) ###
alg_name = "D0CPTagKKPiPi"                       # D0 → K+K−π+π− with a tagged anti-D0
tag_alg = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.773]})   # ψ(3770): 3.773 GeV

# Tag side: the anti-D0, taken from the pre-stored D-tag candidates (one tag_side => single tag)
tag_alg.tag_side(:D0) do |t|
  t.modes :D0toKPiPiPi, :D0toKsPiPi   # standard D-tag modes: Kπππ and KS0 π+π−
  t.charm -1                          # pin the tagged side to the anti-D0
  t.window :deltaE, abs: 0.03         # fully reconstructed tags kept within ~3σ of ΔE (±30 MeV); all other tag-side
                                      # observables (mBC, …) are stored, not cut, and windowed later in ROOT
end

# Signal side: everything the tag did not use — D0 → K+K−π+π−
tag_alg.signal_side do |s|
  s.charged(kp: 1, km: 1, pip: 1, pim: 1)   # exactly one K+, one K−, one π+, one π− (each PID'd under its mass hypothesis)
  s.require_charge 0                        # net charge of the signal side is zero
end

# Kinematic fit: four-momentum conservation + tag D0 mass at its nominal value, χ² < 200
tag_alg.fit do |f|
  f.constrain_four_momentum                                        # 4-momentum conservation against the measured CMS 4-vector
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)     # tag D0 constrained to its nominal mass
  f.chi2_cut 200
end

# Tag-side reconstruction switches (DTagAlg internals stay at their frozen defaults otherwise)
tag_alg.dtag_reconstruction do |d|
  d.species :D0            # D-tag reconstruction for the D0 / anti-D0 species
  d.trim_mode_lists true   # only the declared tag modes are reconstructed
end

# BOSS-side procedures that the tag DSL cannot express, preserved for later use
tag_alg
  .note(:tag_mode_availability,
        "most of the paper's CP-eigenstate tag modes (K+K−, π+π−, K_S0π0π0, K_S0π0, K_S0η, K_S0ω, …) are not present in the
         frozen BOSS D-tag list; only the flavor tags :D0toKPiPiPi (Kπππ) and :D0toKsPiPi (KS0 π+π−) can be declared")
  .note(:tag_side_reconstruction,
        "tag-side sub-detector reconstruction is performed inside the frozen DTagAlg and is not DSL-tunable:
         K_S0 → π+π− with |Vz| < 20 cm, secondary-vertex fit, |M(π+π−) − m(K_S0)| < 12 MeV/c² and decay length > 2x the
         vertex resolution; π0/η from γγ with M(γγ) in [115,150] and [480,580] MeV/c²; η′ via π+π−η (M in [940,976] MeV/c²)
         or ρ0γ (M in [940,970] MeV/c² with M(ππ) in [626,924] MeV/c²); K_L0 modes handled by partial reconstruction;
         Kalman fits constraining the K_S,L0 and D masses for the K_S,L0π+π− tags")
  .note(:background_veto,
        "signal-side π+π− pair required to have an invariant mass outside [477, 507] MeV/c² to veto K_S0 → π+π−
         (window from 3σ of the K_S0 peak); the tag fit builder only offers inside-windows, so this veto is applied at the
         ROOT level on the stored π+π− mass")
  .note(:signal_vertex,
        "signal π+π− pair required to originate from a common vertex within 2x the vertex resolution of the IP;
         no tag-side/signal-side DSL primitive expresses this vertex requirement")

tag_alg.apply   # takes NO Selection argument

# Execute the algorithm on data, inclusive MC and the signal exclusive MC
root_files = tag_alg.execute_on([psi3770_data, psi3770_incMC, exMC_signal])