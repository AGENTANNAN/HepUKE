# Core DSL classes and dependencies are loaded automatically at execution time
### Dataset description ###
# ψ(3770) at 3.773 GeV: 2.93 fb^-1 of real data and the matching inclusive MC
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for Mode I: ψ(3770) → D0 anti-D0, D0 → K- pi+ pi+ pi-
decay_card_modeI = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- pi- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for Mode II: ψ(3770) → D0 anti-D0, D0 → K- pi+ pi0, pi0 → gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples: 200k events for each signal mode
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_D0toK3pi"
  config.related_dataset = psi3770_data
  config.events          = 200000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_D0toKpipi0"
  config.related_dataset = psi3770_data
  config.events          = 200000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ----------------- Mode I: D0 → K- pi+ pi+ pi- -----------------
alg_name_modeI = "D0ToK3pi"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         # K_S0 veto on pi+pi- pairs: flight significance > 2 (no DSL primitive for a flight-significance veto)
         .note(:background_veto, "K_S0 -> pi+pi- veto in Mode I: reject pi+pi- track pairs whose "
                                 "secondary-vertex flight significance (decay length / error) exceeds 2")

sel_modeI = Selection.new
sel_modeI.select_track {   # Charged track selection
            cos_theta 0.93   # |cos(theta)| < 0.93
            Vz        100.0  # |Vz| < 100 cm
            Vr        10.0   # Vr < 10 cm in the transverse plane
            nChrp     ">=2"  # At least two positive tracks
            nChrn     ">=2"  # At least two negative tracks
            nNet      "==0"  # Net charge zero
          }
         .select_photon {  # Photon selection
            tdc_emc_start     0     # TDC start time
            tdc_emc_end       14    # TDC end time
            energyThreshold_b 0.025 # E > 25 MeV in the EMC barrel
            energyThreshold_e 0.050 # E > 50 MeV in the EMC endcap
            angle_to_track    20.0  # Photon > 20 degrees from any charged track
          }
         .pid(method: :probability) {  # Probability-method PID with K/pi separation
            prob_cut 0.001
            identify :kaon, against: [:pion]   # K+ and K- identified against pions
            identify :pion, against: [:kaon]   # pi+ and pi- identified against kaons
            nkm  ">=1"  # At least one K-
            nkp  ">=1"  # At least one K+
            npip ">=2"  # At least two pi+
            npim ">=2"  # At least two pi-
          }
         # Kinematic fit of the K- pi+ pi+ pi- system, mass-constrained to the nominal D0 mass
         .kinematic_fit([:km, :pip, :pip, :pim]) {
            nominal
            invariant_mass_of(:km, :pip, :pip, :pim).constrain_to_nominal_mass_of(:D0)
            chi2_cut 200
          }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ----------------- Mode II: D0 → K- pi+ pi0 -----------------
alg_name_modeII = "D0ToKpiPi0"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {   # Same charged track criteria, relaxed multiplicity
            cos_theta 0.93
            Vz        100.0
            Vr        10.0
            nChrp     ">=1"  # At least one positive track
            nChrn     ">=1"  # At least one negative track
            nNet      "==0"  # Net charge zero
          }
          .select_photon {  # Same photon criteria, at least two photons for the pi0
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track    20.0
            nGam              ">=2"  # At least two photons
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion]
            identify :pion, against: [:kaon]
            nkm  ">=1"
            nkp  ">=1"
            npip ">=1"
            npim ">=1"
          }
          # Kalman fit: reconstruct pi0 from two photons, mass-constrained to the nominal pi0 mass
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"   # At least one pi0 candidate
          }
          # Kinematic fit of the K- pi+ pi0 system, mass-constrained to the nominal D0 mass
          .kinematic_fit([:km, :pip, :pi0]) {
            nominal
            invariant_mass_of(:km, :pip, :pi0).constrain_to_nominal_mass_of(:D0)
            chi2_cut 200
          }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

### Execute on datasets ###
root_files_modeI  = alg_modeI.execute_on([psi3770_data, psi3770_incMC, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([psi3770_data, psi3770_incMC, exMC_modeII])