# BESIII arXiv:1412.4566v2
# Measurement of B(D*0 -> D0 pi0) and B(D*0 -> D0 gamma) at sqrt(s) = 4.009 GeV (482 pb^-1)
# Method: fully reconstruct the D0 anti-D0 pair, mass-constrain both D0 candidates, and
# study the D0 anti-D0 recoil mass to separate D*0 -> D0 pi0 from D*0 -> D0 gamma.
# Five tag combinations (modes I-V) built from three large-BR D0 decay modes.

### Datasets ###
data_4009  = DatasetManager.real_data.find("703_4009")     # 482 pb^-1 at sqrt(s) = 4.009 GeV
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")  # inclusive MC, 500 pb^-1 equivalent

### Decay cards — one per mode; D*0 decays to both D0 pi0 and D0 gamma (PDG fractions) ###
# Mode I: D0 -> K- pi+ , anti-D0 -> K+ pi-
decay_card_I = <<~DECAYCARD
  Decay psi(4260)
  1.000 D*0 anti-D0                               PHSP;
  Enddecay

  Decay D*0
  0.647 D0 pi0                                    PHSP;
  0.353 D0 gamma                                  PHSP;
  Enddecay

  Decay D0
  1.000 K- pi+                                    PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ pi-                                    PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                               PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: D0 -> K- pi+ , anti-D0 -> K+ pi- pi0
decay_card_II = <<~DECAYCARD
  Decay psi(4260)
  1.000 D*0 anti-D0                               PHSP;
  Enddecay

  Decay D*0
  0.647 D0 pi0                                    PHSP;
  0.353 D0 gamma                                  PHSP;
  Enddecay

  Decay D0
  1.000 K- pi+                                    PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ pi- pi0                                PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                               PHSP;
  Enddecay

  End
DECAYCARD

# Mode III: D0 -> K- pi+ pi0 , anti-D0 -> K+ pi-
decay_card_III = <<~DECAYCARD
  Decay psi(4260)
  1.000 D*0 anti-D0                               PHSP;
  Enddecay

  Decay D*0
  0.647 D0 pi0                                    PHSP;
  0.353 D0 gamma                                  PHSP;
  Enddecay

  Decay D0
  1.000 K- pi+ pi0                                PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ pi-                                    PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                               PHSP;
  Enddecay

  End
DECAYCARD

# Mode IV: D0 -> K- pi+ , anti-D0 -> K+ pi- pi+ pi-
decay_card_IV = <<~DECAYCARD
  Decay psi(4260)
  1.000 D*0 anti-D0                               PHSP;
  Enddecay

  Decay D*0
  0.647 D0 pi0                                    PHSP;
  0.353 D0 gamma                                  PHSP;
  Enddecay

  Decay D0
  1.000 K- pi+                                    PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ pi- pi+ pi-                            PHSP;
  Enddecay

  End
DECAYCARD

# Mode V: D0 -> K- pi+ pi+ pi- , anti-D0 -> K+ pi-
decay_card_V = <<~DECAYCARD
  Decay psi(4260)
  1.000 D*0 anti-D0                               PHSP;
  Enddecay

  Decay D*0
  0.647 D0 pi0                                    PHSP;
  0.353 D0 gamma                                  PHSP;
  Enddecay

  Decay D0
  1.000 K- pi+ pi+ pi-                            PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ pi-                                    PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive signal MC — 1,000,000 events in total, split over the five modes ###
exMC_I   = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dst0D0_Kpi_Kpi"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_I
  config.cross_section   = :default
end

exMC_II  = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dst0D0_Kpi_Kpipi0"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_II
  config.cross_section   = :default
end

exMC_III = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dst0D0_Kpipi0_Kpi"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_III
  config.cross_section   = :default
end

exMC_IV  = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dst0D0_Kpi_Kpipipi"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_IV
  config.cross_section   = :default
end

exMC_V   = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dst0D0_Kpipipi_Kpi"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_V
  config.cross_section   = :default
end

### Common notes on the reconstruction method ###
common_notes = <<~NOTES
  The D0 anti-D0 system is the tag; the D*0 decay products (pi0 or gamma) are NOT
  reconstructed. The number of D*0 -> D0 pi0 and D*0 -> D0 gamma events is extracted
  from the square of the D0 anti-D0 recoil mass, whose signal regions are
  [0.01, 0.04] (GeV/c^2)^2 for the pi0 mode and [-0.01, 0.01] (GeV/c^2)^2 for the gamma
  mode. An additional requirement that the momenta of both D0 and anti-D0 be below
  0.65 GeV/c rejects the direct e+e- -> D0 anti-D0 background (peak at 0.75 GeV/c).
  Those requirements act on kinematic-fit-corrected four-momenta and are therefore
  applied in the ROOT analysis stage, not in BOSS.
NOTES

### Algorithm 1: Mode I — D0 -> K- pi+, anti-D0 -> K+ pi- ###
alg_I = Algorithm.new("Dst0D0ModeI")
alg_I.set_header(["Dst0D0ModeIAlg/Dst0D0ModeI.h"])
     .set_constant({"ECMS" => [:double, 4.009]})

sel_I = Selection.new
sel_I.select_track {
      cos_theta 0.93       # |cos(theta)| < 0.93 in the MDC
      Vz        10.0       # within 10 cm of the IP along the beam direction
      Vr        1.0        # within 1 cm of the IP in the transverse plane
      nChrp     "==2"      # K+ and pi+ (one oppositely charged kaon/pion pair each side)
      nChrn     "==2"      # K- and pi-
      nNet      "==0"
    }
    .select_photon {
      angle_to_track    20.0   # photon at least 20 degrees away from any charged track
      energyThreshold_b 0.025  # > 25 MeV in the EMC barrel
      energyThreshold_e 0.050  # > 50 MeV in the EMC end-cap
      tdc_emc_start     0      # EMC time 0-700 ns, coincident with the collision
      tdc_emc_end       14
    }
    .pid(method: :probability) {
      prob_cut 0.001   # P_pi (P_K) > 0.1% and P_pi > P_K (P_K > P_pi) from TOF and dE/dx
      identify :kaon, against: [:pion]
      identify :pion, against: [:kaon]
      nkp  "==1"
      nkm  "==1"
      npip "==1"
      npim "==1"
    }
    # Fit with the D0 and anti-D0 candidates constrained to the nominal D0 mass
    .kinematic_fit([:km, :pip, :kp, :pim]) {
      nominal
      invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)
      invariant_mass_of(:kp, :pim).constrain_to_nominal_mass_of(:D0)
      chi2_cut 30
    }

alg_I.note(:method,
  "Reconstruct the D0 anti-D0 pair and use the mass recoiling against it to separate " \
  "D*0 -> D0 pi0 from D*0 -> D0 gamma. The pi0 / gamma from the D*0 decay is not " \
  "reconstructed (partial reconstruction of the D*0).")
     .note(:recoil_mass_regions,
  common_notes)
     .note(:multi_candidate,
  "When more than one D0 anti-D0 combination satisfies the selection criteria, the " \
  "combination with the smallest total chi2 is selected.")
     .note(:helix_correction,
  "Corrected track parameters are used in the nominal MC simulation; the difference in " \
  "the measured branching fractions with and without this correction is taken as the " \
  "systematic uncertainty of the kinematic-fit chi2 requirement.")
     .note(:fsr_simulation,
  "The fraction of events with final-state-radiation photons from charged pions in data " \
  "is 20% higher than in MC; the FSR fraction in MC is enlarged by 1.2^X (X = number of " \
  "charged pions) to estimate the associated systematic uncertainty.")

alg_I.with_decay_card(decay_card_I).apply(sel_I)
alg_I.execute_on([data_4009, incMC_4009, exMC_I])

### Algorithm 2: Mode II — D0 -> K- pi+, anti-D0 -> K+ pi- pi0 ###
alg_II = Algorithm.new("Dst0D0ModeII")
alg_II.set_header(["Dst0D0ModeIIAlg/Dst0D0ModeII.h"])
      .set_constant({"ECMS" => [:double, 4.009]})

sel_II = Selection.new
sel_II.select_track {
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==2"     # K+ pi+ (anti-D0 side pi- is also negative, so 2 pi+ / 2 pi-)
       nChrn     "==2"
       nNet      "==0"
     }
     .select_photon {
       angle_to_track    20.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       tdc_emc_start     0
       tdc_emc_end       14
       nGam              ">=2"   # at least two good photons for the pi0
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon, against: [:pion]
       identify :pion, against: [:kaon]
       nkp  "==1"
       nkm  "==1"
       npip "==1"
       npim "==1"
     }
     # pi0 from the anti-D0 -> K+ pi- pi0 decay, mass-constrained
     .kalman_kinematic_fit([:gamma, :gamma]) {
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 25
       npi0 ">=1"
     }
     # Fit constraining the D0 (K- pi+) and the anti-D0 (K+ pi- pi0) to the nominal D0 mass
     .kinematic_fit([:km, :pip, :kp, :pim, :pi0]) {
       nominal
       invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)
       invariant_mass_of(:kp, :pim, :pi0).constrain_to_nominal_mass_of(:D0)
       chi2_cut 30
     }

alg_II.note(:pi0_preselection,
  "The invariant mass of the two photons is required to be within +/-15 MeV/c^2 of the " \
  "nominal pi0 mass before the pi0 mass constraint is included in the kinematic fit.")
      .note(:method,
  "D0 anti-D0 recoil-mass method; the pi0 from the D*0 decay is not reconstructed.")
      .note(:recoil_mass_regions,
  common_notes)
      .note(:multi_candidate,
  "Smallest total chi2 combination retained.")
      .note(:helix_correction,
  "Track helix-parameter correction applied in the nominal MC (as in the other modes).")

alg_II.with_decay_card(decay_card_II).apply(sel_II)
alg_II.execute_on([data_4009, incMC_4009, exMC_II])

### Algorithm 3: Mode III — D0 -> K- pi+ pi0, anti-D0 -> K+ pi- ###
alg_III = Algorithm.new("Dst0D0ModeIII")
alg_III.set_header(["Dst0D0ModeIIIAlg/Dst0D0ModeIII.h"])
       .set_constant({"ECMS" => [:double, 4.009]})

sel_III = Selection.new
sel_III.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
      }
      .select_photon {
        angle_to_track    20.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        tdc_emc_start     0
        tdc_emc_end       14
        nGam              ">=2"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion]
        identify :pion, against: [:kaon]
        nkp  "==1"
        nkm  "==1"
        npip "==1"
        npim "==1"
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      # D0 (K- pi+ pi0) and anti-D0 (K+ pi-) both constrained to the nominal D0 mass
      .kinematic_fit([:km, :pip, :pi0, :kp, :pim]) {
        nominal
        invariant_mass_of(:km, :pip, :pi0).constrain_to_nominal_mass_of(:D0)
        invariant_mass_of(:kp, :pim).constrain_to_nominal_mass_of(:D0)
        chi2_cut 30
      }

alg_III.note(:pi0_preselection,
  "The two-photon invariant mass must lie within +/-15 MeV/c^2 of the nominal pi0 mass.")
       .note(:method,
  "D0 anti-D0 recoil-mass method; the pi0 from the D*0 decay is not reconstructed.")
       .note(:recoil_mass_regions,
  common_notes)
       .note(:multi_candidate,
  "Smallest total chi2 combination retained.")

alg_III.with_decay_card(decay_card_III).apply(sel_III)
alg_III.execute_on([data_4009, incMC_4009, exMC_III])

### Algorithm 4: Mode IV — D0 -> K- pi+, anti-D0 -> K+ pi- pi+ pi- ###
alg_IV = Algorithm.new("Dst0D0ModeIV")
alg_IV.set_header(["Dst0D0ModeIVAlg/Dst0D0ModeIV.h"])
      .set_constant({"ECMS" => [:double, 4.009]})

sel_IV = Selection.new
sel_IV.select_track {
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==3"     # K+ pi+ pi+ from the anti-D0 four-pion side plus pi+ from K- pi+
       nChrn     "==3"     # K- pi- pi-
       nNet      "==0"
     }
     .select_photon {
       angle_to_track    20.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       tdc_emc_start     0
       tdc_emc_end       14
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon, against: [:pion]
       identify :pion, against: [:kaon]
       nkp  "==1"
       nkm  "==1"
       npip "==2"   # two oppositely charged pion pairs in the final state
       npim "==2"
     }
     .kinematic_fit([:km, :pip, :kp, :pim, :pip, :pim]) {
       nominal
       invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)
       invariant_mass_of(:kp, :pim, :pip, :pim).constrain_to_nominal_mass_of(:D0)
       chi2_cut 30
     }

alg_IV.note(:method,
  "Modes IV and V contain one oppositely charged kaon pair and two oppositely charged " \
  "pion pairs; combinations with more than six charged tracks are not used.")
       .note(:recoil_mass_regions,
  common_notes)
       .note(:multi_candidate,
  "Smallest total chi2 combination retained.")

alg_IV.with_decay_card(decay_card_IV).apply(sel_IV)
alg_IV.execute_on([data_4009, incMC_4009, exMC_IV])

### Algorithm 5: Mode V — D0 -> K- pi+ pi+ pi-, anti-D0 -> K+ pi- ###
alg_V = Algorithm.new("Dst0D0ModeV")
alg_V.set_header(["Dst0D0ModeVAlg/Dst0D0ModeV.h"])
     .set_constant({"ECMS" => [:double, 4.009]})

sel_V = Selection.new
sel_V.select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==3"
      nChrn     "==3"
      nNet      "==0"
    }
    .select_photon {
      angle_to_track    20.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      tdc_emc_start     0
      tdc_emc_end       14
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion]
      identify :pion, against: [:kaon]
      nkp  "==1"
      nkm  "==1"
      npip "==2"
      npim "==2"
    }
    .kinematic_fit([:km, :pip, :pip, :pim, :kp, :pim]) {
      nominal
      invariant_mass_of(:km, :pip, :pip, :pim).constrain_to_nominal_mass_of(:D0)
      invariant_mass_of(:kp, :pim).constrain_to_nominal_mass_of(:D0)
      chi2_cut 30
    }

alg_V.note(:method,
  "D0 anti-D0 recoil-mass method; the pi0 or gamma from the D*0 decay is not reconstructed.")
     .note(:recoil_mass_regions,
  common_notes)
     .note(:multi_candidate,
  "Smallest total chi2 combination retained.")
     .note(:background,
  "The dominant backgrounds are open-charm processes (cross section 7.1 nb) and ISR " \
  "production of psi(3770) with psi(3770) -> D0 anti-D0 (cross section 0.114 nb); their " \
  "numbers are taken from the inclusive MC sample.")

alg_V.with_decay_card(decay_card_V).apply(sel_V)
alg_V.execute_on([data_4009, incMC_4009, exMC_V])
