### Dataset preparation ###
# psi(3770) at 3.773 GeV : real data (7.93 fb^-1) + inclusive MC
d3773_data  = DatasetManager.real_data.find("712_3773")
d3773_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Signal decay card: psi(3770) -> D+ D-, D+ -> K- pi+ pi+ pi0 (pi0 -> gamma gamma), tag D- -> K+ pi- pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D+  D-                PHSP;
    Enddecay

    Decay D+
    1.0000  K-  pi+  pi+  pi0     PHSP;
    Enddecay

    Decay D-
    1.0000  K+  pi-  pi-          PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma          PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC generated flat in phase space (reweighted later by the amplitude-analysis result)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dp_Kpipipi0_Dm_Kpipi"
  config.related_dataset = d3773_data       # associated real dataset
  config.events          = 1_000_000        # 1e6 events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "DstagDpKpipipi0"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })
   .set_alias({ "std::vector<double>" => "Vdouble" })

# Tag side: D- -> K+ pi- pi-  (D+- family, charm pinned to -1 for the D- tag).
# The tag candidate carries its own selected and PID'd tracks/showers.
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi
  t.charm -1
  # Explicit tag-side DeltaE window requested by the analysis.
  t.window :deltaE, min: -0.025, max: 0.024
end

# Signal side: fully reconstructed D+ -> K- pi+ pi+ pi0 (pi0 -> gamma gamma),
# built from the tracks and showers the tag did not use.
alg.signal_side do |s|
  s.charged(km: 1, pip: 2)      # one K- and two pi+ ...
  s.require_charge 1            # ... with total charge +1
  s.photons 2                   # two photons from pi0
  s.min_photon_angle 10.0       # opening angle to charged tracks > 10 deg
  s.min_photon_energy 0.025     # barrel shower energy floor (> 25 MeV)
end

# 6C kinematic fit: four-momentum conservation (4C) + pi0 mass (1C) + D tag mass (1C)
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)                 # pi0 mass window
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)    # pi0 mass constraint
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)           # D tag mass constraint
  f.chi2_cut 100
end

# BOSS-side procedures that cannot be expressed in the tag DSL
alg.note(:pi0_reconstruction,
         "pi0 reconstructed from gamma gamma through a separate 1C mass-constrained Kalman fit to the
          nominal pi0 mass with chi2 < 30, applied before the 6C fit; the gamma gamma invariant mass
          is required in [0.115, 0.150] GeV/c^2 and at least one photon must lie in the barrel EMC")
alg.note(:efficiency_curve,
         "exclusive signal MC is generated flat in phase space; it is to be reweighted by the
          amplitude-analysis result to obtain the reconstruction efficiency")

# Validate / render the tag spec (no Selection argument), then run
alg.apply
root_files = alg.execute_on([d3773_data, d3773_incMC, exMC_signal])