# ============================================================
#  Datasets — signal resonance, psi(3770) and off-resonance points
# ============================================================
# psi(2S) at 3.686 GeV (signal)
psip_data     = DatasetManager.real_data.find("709_3686")
psip_incMC    = DatasetManager.inclusive_mc.find("709_3686")

# psi(3770) at 3.773 GeV
psipp_data    = DatasetManager.real_data.find("712_3773")
psipp_incMC   = DatasetManager.inclusive_mc.find("712_3773")

# Off-resonance continuum points
off3650_data  = DatasetManager.real_data.find("709_3650")   # 3.650 GeV
off3650_incMC = DatasetManager.inclusive_mc.find("709_3650")
off3682_data  = DatasetManager.real_data.find("709_3682")   # 3.682 GeV
off3682_incMC = DatasetManager.inclusive_mc.find("709_3682")

# ============================================================
#  Signal decay card: e+e- -> gamma pi0, pi0 -> gamma gamma
#  (psi(4260) is the BESIII/KKMC top-mother convention)
# ============================================================
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for e+e- -> gamma pi0
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psi2S_gammapi0"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# ============================================================
#  Algorithm
# ============================================================
alg_name  = "GammaPi0"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

# ============================================================
#  Event selection (BOSS) — three-photon final state
# ============================================================
event_selection = Selection.new
event_selection
  # Veto all charged tracks: no good track may survive the quality cuts
  .select_track {
    cos_theta 0.93    # |cos(theta)| < 0.93 for charged tracks
    Vz        10.0    # |Vz| < 10 cm
    Vr        1.0     # Vr < 1 cm
    nTot      "==0"   # exactly zero good charged tracks
  }
  # Exactly three photons with the standard timing / shower-quality cuts
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0    # shower must be > 10 deg from any charged track
    energyThreshold_b 0.025   # > 25 MeV in the EMC barrel (|cos(theta_gamma)| < 0.8)
    energyThreshold_e 0.050
    nGam              "==3"
  }
  # 4C kinematic fit over the three photons (no pi0 mass constraint)
  .kinematic_fit([:gamma, :gamma, :gamma]) {
    nominal                   # nominal fit — corrected four-momenta are saved
    constrain_four_momentum   # 4C energy-momentum conservation to the CMS energy
    chi2_cut 40               # chi2 < 40
  }

# BOSS-side procedures that have no dedicated DSL construct
algorithm
  .note(:photon_selection,
        "all three photons are required to lie in the EMC barrel, " \
        "|cos(theta_gamma)| < 0.8, in addition to the standard shower-quality " \
        "and timing cuts; no photon isolation cut is applied")
  .note(:photon_conversion_veto,
        "each photon must have fewer than 8 MDC hits along the straight line " \
        "from the interaction point to the matched EMC shower, to reject " \
        "Converted photons")
  .note(:background_veto,
        "the highest-energy photon is identified as the radiative photon and " \
        "the remaining two as the pi0 daughters; the helicity angle " \
        "cos(theta_hel) < 0.7 suppresses the e+e- -> gamma gamma gamma QED " \
        "background")
  .note(:beam_energy,
        "the datasets span several CMS energies (3.686, 3.773, 3.650, 3.682 " \
        "GeV); ECMS is set to the psi(2S) nominal value and the per-dataset " \
        "beam energy must be used by the 4C kinematic fit")

# Attach decay card and render the selection into the BOSS algorithm
algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# ============================================================
#  Execute on real data, inclusive MC and signal exclusive MC
# ============================================================
root_files = algorithm.execute_on([
  psip_data,   psip_incMC,
  psipp_data,  psipp_incMC,
  off3650_data, off3650_incMC,
  off3682_data, off3682_incMC,
  exMC_signal
])