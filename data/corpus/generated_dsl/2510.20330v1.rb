# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
# √s = 4.178 GeV, 3.19 fb⁻¹ → BOSS 703, sample 703_4180 (4178 MeV, L = 3189.0 pb⁻¹)
data_4178  = DatasetManager.real_data.find("703_4180")       # Real data at 4.178 GeV
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")    # Corresponding inclusive MC

# Decay card for the signal chain: e+e- -> Ds*+ Ds-, Ds*+ -> Ds+ pi0, Ds+ -> K+ K- pi+
# (charge conjugate implied). psi(4260) is the KKMC top-mother convention.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 Ds*+ Ds- PHSP;
    Enddecay

    Decay Ds*+
    1.0 Ds+ pi0 PHSP;
    Enddecay

    Decay Ds+
    1.0 K+ K- pi+ PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the calibration chain: e+e- -> D*+ D-, D*+ -> D+ pi0, D+ -> K- pi+ pi+
decay_card_calib = <<~DECAYCARD
    Decay psi(4260)
    1.0 D*+ D- PHSP;
    Enddecay

    Decay D*+
    1.0 D+ pi0 PHSP;
    Enddecay

    Decay D+
    1.0 K- pi+ pi+ PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for each of the two decay modes
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_dsstar_dspi0_kkpi"
  config.related_dataset = data_4178
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exMC_calib = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_dstar_dpi0_kpipi"
  config.related_dataset = data_4178
  config.events          = 100000
  config.decay_card      = decay_card_calib
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ---------------------------------------------------------------------------
# Signal channel: Ds*+ -> Ds+ pi0, Ds+ -> K+ K- pi+   (two independent modes →
# separate Algorithm objects per Rule T1)
# ---------------------------------------------------------------------------
alg_name_sig = "DsStarDsPi0"
alg_sig = Algorithm.new(alg_name_sig)
alg_sig.set_header(["#{alg_name_sig}Alg/#{alg_name_sig}.h"])
       .set_constant({"ECMS" => [:double, 4.178]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_sig = Selection.new
sel_sig.select_track {                     # Charged track selection
           cos_theta 0.93                  # |cos(theta)| < 0.93
           Vz        10.0                  # |Vz| < 10 cm
           Vr        1.0                   # Vr < 1 cm
           nChrp     ">=2"                 # at least two positive tracks
           nChrn     ">=1"                 # at least one negative track
         }
       .select_photon {                    # Photon selection
           tdc_emc_start     0             # EMC TDC window start
           tdc_emc_end       14            # EMC TDC window end
           energyThreshold_b 0.025         # > 25 MeV in the barrel
           energyThreshold_e 0.050         # > 50 MeV in the endcap
           nGam              ">=2"         # at least two photons
         }
       .pid(method: :probability) {        # Probability-based PID with K/pi separation
           prob_cut 0.001                  # selection probability > 0.001
           identify :kaon, against: [:pion, :proton]   # K+ / K- candidates
           identify :pion, against: [:kaon, :proton]   # pi+ / pi- candidates
           nkp  ">=1"                      # at least one K+
           nkm  ">=1"                      # at least one K-
           npip ">=1"                      # at least one pi+
         }
       .kalman_kinematic_fit([:gamma, :gamma]) {       # Reconstruct pi0 from gamma gamma
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 25                     # chi2 < 25
           npi0 ">=1"                      # at least one pi0
         }
       .kinematic_fit([:kp, :km, :pip, :pi0]) {        # 4C fit on K+ K- pi+ pi0
           nominal
           constrain_four_momentum         # constrain total 4-momentum to the CMS
           chi2_cut 200                    # loose BOSS chi2 cut (tight cut applied in ROOT)
         }

alg_sig.with_decay_card(decay_card_signal).apply(sel_sig)

# ---------------------------------------------------------------------------
# Calibration channel: D*+ -> D+ pi0, D+ -> K- pi+ pi+
# ---------------------------------------------------------------------------
alg_name_cal = "DStarDPi0"
alg_cal = Algorithm.new(alg_name_cal)
alg_cal.set_header(["#{alg_name_cal}Alg/#{alg_name_cal}.h"])
       .set_constant({"ECMS" => [:double, 4.178]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_cal = Selection.new
sel_cal.select_track {                     # Same charged track quality cuts
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp     ">=2"                 # at least two positive tracks
           nChrn     ">=1"                 # at least one negative track
         }
       .select_photon {                    # Same photon selection
           tdc_emc_start     0
           tdc_emc_end       14
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam              ">=2"
         }
       .pid(method: :probability) {        # Same probability-based PID
           prob_cut 0.001
           identify :kaon, against: [:pion, :proton]
           identify :pion, against: [:kaon, :proton]
           nkm  ">=1"                      # at least one K-
           npip ">=2"                      # at least two pi+
         }
       .kalman_kinematic_fit([:gamma, :gamma]) {       # Reconstruct pi0 from gamma gamma
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 25
           npi0 ">=1"
         }
       .kinematic_fit([:km, :pip, :pip, :pi0]) {       # 4C fit on K- pi+ pi+ pi0
           nominal
           constrain_four_momentum
           chi2_cut 200
         }

alg_cal.with_decay_card(decay_card_calib).apply(sel_cal)

# Execute the two algorithms on real data, inclusive MC, and the corresponding exclusive MC
root_files_sig = alg_sig.execute_on([data_4178, incMC_4178, exMC_signal])
root_files_cal = alg_cal.execute_on([data_4178, incMC_4178, exMC_calib])