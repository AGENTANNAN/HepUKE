# Core DSL classes and dependencies will be loaded automatically at execution.

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # tag-based reformed inclusive MC at 3.773 GeV

# Decay card for the signal process: psi(3770) -> D0 D0bar (quantum-correlated pair)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 1,000,000 psi(3770) -> D0 D0bar events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0D0bar"
  config.related_dataset = data_3773          # associated real dataset
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
# Quantum-correlated D0 D0bar double-tag analysis (tag1 = flavor/CP tags, tag2 = signal mode)
alg_name = "D0D0barKsPiPiDT"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})   # 3.773 GeV centre-of-mass energy
   .note(:tag_mode_availability, "Semileptonic tags together with the K_L0 pi0, K_L0 pi0 pi0 and
     several CP-odd tag modes are not available in the DTagAlg tag collection; their single-tag
     yields are therefore taken from the known branching fractions instead of being reconstructed.")
   .note(:kl_pipi_missing_mass, "The K_L0 pi+ pi- signal is not accessed through a reconstructed
     tag but is reconstructed via the missing-mass-squared variable.")
   .with_decay_card(decay_card_signal)           # variables follow from psi(3770) -> D0 D0bar

# --- Tag side 1: flavor, CP-even, CP-odd and mixed-CP tags (both tag charges scanned) ---
alg.tag_side(:D0) do |t|
  t.modes(:D0toKPi,       # K-+ pi+-
          :D0toKPiPi0,    # K-+ pi+- pi0
          :D0toKPiPiPi,   # K-+ pi+- pi+-
          :D0toKK,        # K+ K-
          :D0toPiPi,      # pi+ pi-
          :D0toKsPi0,     # K_S0 pi0
          :D0toPiPiPi0,   # pi+ pi- pi0
          :D0toKsEta,     # K_S0 eta
          :D0toKsOmega,   # K_S0 omega
          :D0toKsPiPi)    # K_S0 pi+ pi-
end

# --- Tag side 2: the signal mode D0 -> K_S0 pi+ pi- (two tag_side calls => double tag) ---
alg.tag_side(:D0) do |t|
  t.modes(:D0toKsPiPi)
end

# --- Signal side: everything is consumed by the two tags; no additional photons or tracks ---
alg.signal_side do |s|
  s.photons 0            # no extra photons on the (empty) non-tag side
                         # charged omitted -> zero leftover charged tracks required
end

# --- Kinematic fit: four-momentum conservation + both tag masses at the nominal D0 mass ---
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# Render the tag analysis (no Selection argument) and run on data, inclusive MC and signal MC
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_signal])