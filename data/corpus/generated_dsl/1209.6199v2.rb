# Core DSL classes and dependencies are loaded automatically at execution.
### Dataset description ###
# ψ(2S) on-resonance data collected in 2009 and its matching inclusive MC
psip_data    = DatasetManager.real_data.find("709_3686")
psip_incMC   = DatasetManager.inclusive_mc.find("709_3686")
# 3.650 GeV off-resonance data (~44 pb^-1) and its continuum inclusive MC
offres_data  = DatasetManager.real_data.find("709_3650")
offres_incMC = DatasetManager.inclusive_mc.find("709_3650")

# Exclusive MC for the inclusive hadronic psi(2S) decay, generated with PHSP.
# A set of representative hadronic final states stands in for the inclusive
# hadron sample used to estimate the detection efficiency.
decay_card_hadrons = <<~DECAYCARD
    Decay psi(2S)
    0.25 pi+ pi-             PHSP;
    0.25 pi+ pi- pi0         PHSP;
    0.25 pi+ pi- pi+ pi-     PHSP;
    0.25 K+ K- pi0           PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 10^6-event exclusive psi(2S) -> anything (phase-space) sample for efficiency
exMC_hadrons = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_inclusive_hadrons"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_hadrons
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipInclusiveHadrons"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})   # psi(2S) center-of-mass energy in GeV

event_selection = Selection.new
event_selection
  .select_track {                       # Charged-track quality selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        15.0                      # |Vz| < 15 cm
    Vr        1.0                       # Vr < 1 cm
    nTot      ">=1"                     # at least one good track
  }
  .select_photon {                      # Photon (EMC shower) selection
    tdc_emc_start     0                 # TDC start time
    tdc_emc_end       14                # TDC end time
    angle_to_track    10.0              # separated from any track by > 10 degrees
    energyThreshold_b 0.025             # 25 MeV in the barrel (|cos(theta)| < 0.8)
    energyThreshold_e 0.050             # 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
    nGam              ">=2"             # at least two photons
  }

# No PID and no kinematic fit are applied in this inclusive-hadron counting analysis.

# BOSS-side steps that have no dedicated DSL construct are preserved as notes.
alg
  .note(:background_veto,
        "For events with exactly two good charged tracks (dominated by Bhabha and dimuon
         events) each track is required to have momentum < 1.7 GeV/c and the opening
         angle between the two tracks < 176 deg; events with more than two good tracks
         pass unchanged as the inclusive-hadron topology.")
  .note(:one_track_pi0_requirement,
        "For events with exactly one good charged track, at least two additional photons
         are required and the photon pair whose gamma-gamma invariant mass is closest to
         the pi0 is selected, keeping only |M(gamma gamma) - M(pi0)| < 0.015 GeV/c^2.")
  .note(:visible_energy_cut,
        "Event-level visible-energy cut E_vis / E_cm > 0.4, applied after track and
         photon selection.")
  .note(:vertex_z_signal_sideband,
        "Average per-track z vertex |Vbar_z| is used to separate signal
         (|Vbar_z| < 4.0 cm) from a non-collision-background sideband
         (6.0 < |Vbar_z| < 10.0 cm); N_obs = N_signal - N_sideband, followed by
         subtraction of the off-resonance continuum with scaling factor f = 3.677.")

# Render the algorithm and the selection chain for the given process
alg.with_decay_card(decay_card_hadrons).apply(event_selection)

# Execute on real data, inclusive MC and the exclusive efficiency sample
root_files = alg.execute_on([psip_data, psip_incMC, offres_data, offres_incMC, exMC_hadrons])