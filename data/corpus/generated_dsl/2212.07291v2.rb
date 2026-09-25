# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description (nine CMS energy points) ###
# 4.661, 4.682, 4.699, 4.740, 4.750, 4.781, 4.843, 4.918, 4.951 GeV
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")

incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_4740 = DatasetManager.inclusive_mc.find("707_4740")
incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

energy_points  = [data_4660, data_4680, data_4700, data_4740, data_4750,
                  data_4780, data_4840, data_4914, data_4946]
incMC_samples  = [incMC_4660, incMC_4680, incMC_4700, incMC_4740, incMC_4750,
                  incMC_4780, incMC_4840, incMC_4914, incMC_4946]

# Decay card: e+e- -> omega X(3872), X(3872) -> pi+pi- J/psi, J/psi -> e+e-,
#             omega -> pi+pi-pi0, pi0 -> gamma gamma (top mother = psi(4260), KKMC convention)
decay_card_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: same signal but J/psi -> mu+mu-
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive MC per energy point, for each lepton channel (e and mu)
exMCs_ee = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_omegaX3872_ee"
  config.events        = 200_000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

exMCs_mumu = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_omegaX3872_mumu"
  config.events        = 200_000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ---------------------------------------------------------------------------
# Topology I — 6-track (all-charged) chain: the pi0 is not reconstructed,
# its recoil mass is constrained in a 1C fit.
# ---------------------------------------------------------------------------
alg_name_6trk = "OmegaX3872SixTrk"
alg_6trk = Algorithm.new(alg_name_6trk)
alg_6trk.set_header(["#{alg_name_6trk}Alg/#{alg_name_6trk}.h"])
        .set_constant({"ECMS" => [:double, 4.661]})

sel_6trk = Selection.new
sel_6trk.select_track {                 # charged tracks: 6 tracks, net charge 0
            cos_theta 0.93              # |cos(theta)| < 0.93
            Vz        10.0              # |Vz| < 10 cm
            Vr        1.0               # Vr < 1 cm
            nChrp     "==3"
            nChrn     "==3"
            nNet      "==0"
        }
        .select_photon {                # photons: TDC 0-14, 25/50 MeV, >10 deg from tracks
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track    10.0
            nGam              ">=1"
        }
        .pid(method: :probability) {    # probability PID: leptons vs pi/K, pions vs K
            prob_cut 0.001
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 0.6
            identify :pion, against: [:kaon]
            nlp "==1"                   # the J/psi lepton pair (l+ l-)
            nlm "==1"
        }
        # 1C kinematic fit: recoil mass of 4pi + ll against e+e- constrained to pi0 mass
        .kinematic_fit([:pip, :pim, :pip, :pim, :lp, :lm]) {
            nominal
            miss_track_of :pi0          # undetected pi0 (recoil constrained to pi0 mass)
            chi2_cut 15
        }

# ---------------------------------------------------------------------------
# Topology II — 5-track chain: one pion is missed; pi0 reconstructed from gamma gamma.
# ---------------------------------------------------------------------------
alg_name_5trk = "OmegaX3872FiveTrkMissPi"
alg_5trk = Algorithm.new(alg_name_5trk)
alg_5trk.set_header(["#{alg_name_5trk}Alg/#{alg_name_5trk}.h"])
        .set_constant({"ECMS" => [:double, 4.661]})

sel_5trk = Selection.new
sel_5trk.select_track {                 # charged tracks: 5 tracks
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nTot      "==5"
        }
        .select_photon {                # photons for the pi0 -> gamma gamma
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track    10.0
            nGam              ">=2"
        }
        .pid(method: :probability) {    # probability PID: leptons vs pi/K, pions vs K
            prob_cut 0.001
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 0.6
            identify :pion, against: [:kaon]
            nlp "==1"
            nlm "==1"
        }
        # pi0 reconstruction from gamma gamma with a mass constraint
        .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 200
            npi0 ">=1"
        }
        # 2C kinematic fit: gamma gamma -> pi0 mass + recoil mass of 3pi + pi0 + ll to pi mass
        .kinematic_fit([:pip, :pim, :pip, :gamma, :gamma, :lp, :lm]) {
            nominal
            miss_track_of :pim          # the missed pion (recoil constrained to pi mass)
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
        }

# Notes on BOSS-side procedures that are not fully expressible in the DSL
[alg_6trk, alg_5trk].each do |alg|
  alg.note(:pid_correction_method,
           "J/psi lepton pair chosen as the two highest-momentum (>1 GeV/c) opposite-charge " \
           "tracks; electrons identified by EMC energy > 0.8 GeV, muons by EMC energy < 0.4 GeV " \
           "requiring at least one MUC layer with penetration depth > 3. Not fully expressible " \
           "via identify_high_momentum_leptons (which uses fixed v1 thresholds).")
     .note(:energy_dependent_ecms,
           "Analysis spans nine CMS energy points (4.661, 4.682, 4.699, 4.740, 4.750, 4.781, " \
           "4.843, 4.918, 4.951 GeV); ECMS is set to the first point and must be reset per " \
           "energy when running the scan.")
end

# Attach the decay card defining the kinematic variables, and apply the selection
alg_6trk.with_decay_card(decay_card_ee).apply(sel_6trk)
alg_5trk.with_decay_card(decay_card_ee).apply(sel_5trk)

# Execute on all real data, inclusive MC, and signal MC (both lepton channels)
all_datasets = energy_points + incMC_samples + exMCs_ee + exMCs_mumu
root_files_6trk = alg_6trk.execute_on(all_datasets)
root_files_5trk = alg_5trk.execute_on(all_datasets)