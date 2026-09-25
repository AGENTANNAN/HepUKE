# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# e+e- -> K+ K- J/psi (J/psi -> l+ l-, l = e or mu) at twelve CMS energies,
# 4.612 -- 4.951 GeV (BOSS releases 706 and 707).
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")

incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_4740 = DatasetManager.inclusive_mc.find("707_4740")
incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

data_points  = [data_4610, data_4620, data_4640, data_4660, data_4680, data_4700,
                data_4740, data_4750, data_4780, data_4840, data_4914, data_4946]
incmc_points = [incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700,
                incMC_4740, incMC_4750, incMC_4780, incMC_4840, incMC_4914, incMC_4946]

# Decay card for the signal process e+e- -> K+ K- J/psi, J/psi -> mu+ mu-
# (the same signal MC serves both the full- and partial-reconstruction channels).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000  K+  K-  J/psi     PHSP;
    Enddecay

    Decay J/psi
    1.0000  mu+  mu-         PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive signal MC at every energy point (one sample per point)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_KKJpsi_mumu"   # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ---------------- Full reconstruction: e+e- -> K+ K- l+ l- ----------------
alg_full_name = "KKJpsiFull"
alg_full = Algorithm.new(alg_full_name)
alg_full.set_header(["#{alg_full_name}Alg/#{alg_full_name}.h"])
        .set_constant({"ECMS" => [:double, 4.78]})  # representative; per-point energy set at execution

full_selection = Selection.new
full_selection
    .select_track {                     # charged tracks
        cos_theta 0.93                  # |cos(theta)| < 0.93
        Vz        10.0                  # |Vz| < 10 cm
        Vr        1.0                   # Vr < 1 cm
        nChrp     "==2"                 # K+ and l+
        nChrn     "==2"                 # K- and l-
        nNet      "==0"
    }
    .pid(method: :probability) {        # lepton + kaon identification
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.95,
                                       treat_as_electron_if_energy_above: 0.6
        identify :kaon, against: [:pion]   # K+ and K- against pions
        nkp "==1"
        nkm "==1"
        nlp "==1"                       # 1 l+
        nlm "==1"                       # 1 l-
    }
    .kinematic_fit([:kp, :km, :lp, :lm]) {   # 4C fit to K+ K- l+ l-
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

alg_full
    .note(:energy_dependent_pid_cut, "high-momentum lepton threshold is p > 0.95 GeV/c for sqrt(s) < 4.84 GeV and p > 1.05 GeV/c for sqrt(s) >= 4.84 GeV; the DSL high-momentum lepton cut uses a single threshold, retune per energy point")
    .note(:emu_separation, "electron/muon separation via EMC energy applied after identification: muon E_EMC < 0.4 GeV, electron E_EMC > 1.0 GeV")
    .note(:muon_muc_depth_cut, "mu+mu- mode: require at least one muon with MUC penetration depth > 30 cm")
    .note(:background_veto, "e+e- mode: require cos(theta)(K+K-) < 0.98 to reject radiative Bhabha events")
    .note(:ecms_scan, "ECMS varies over the twelve scan points; the kinematic fit must use the per-point collision energy")
    .with_decay_card(decay_card_signal)
    .apply(full_selection)

# ---------------- Partial reconstruction: e+e- -> K+ l+ l- (K- missing) ----------------
alg_partial_name = "KKJpsiPartial"
alg_partial = Algorithm.new(alg_partial_name)
alg_partial.set_header(["#{alg_partial_name}Alg/#{alg_partial_name}.h"])
           .set_constant({"ECMS" => [:double, 4.78]})  # representative; per-point energy set at execution

partial_selection = Selection.new
partial_selection
    .select_track {                     # charged tracks (K- is not reconstructed)
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"                 # K+ and l+
        nChrn     "==1"                 # l-
        nNet      "==1"
    }
    .pid(method: :probability) {        # lepton + kaon identification
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.95,
                                       treat_as_electron_if_energy_above: 0.6
        identify :kaon, against: [:pion]   # charged kaon against pions
        nkp "==1"
        nlp "==1"                       # 1 l+
        nlm "==1"                       # 1 l-
    }
    .kinematic_fit([:kp, :lp, :lm]) {   # 1C fit to K+ l+ l- with the K- missing
        nominal
        miss_track_of :km
        constrain_four_momentum
        chi2_cut 200
    }

alg_partial
    .note(:energy_dependent_pid_cut, "high-momentum lepton threshold is p > 0.95 GeV/c for sqrt(s) < 4.84 GeV and p > 1.05 GeV/c for sqrt(s) >= 4.84 GeV; the DSL high-momentum lepton cut uses a single threshold, retune per energy point")
    .note(:emu_separation, "electron/muon separation via EMC energy applied after identification: muon E_EMC < 0.4 GeV, electron E_EMC > 1.0 GeV")
    .note(:muon_muc_depth_cut, "mu+mu- mode: require both muons to have MUC penetration depth > 30 cm")
    .note(:ee_mode_angular_cuts, "e+e- mode: cos(theta)(e+) < 0.8, cos(theta)(e-) > -0.8, |cos(theta)(K+-)| < 0.8, |cos(alpha_K K_miss)| < 0.95, |cos(alpha_K e+-)| < 0.95 and |cos(alpha_K e-+)| < 0.95")
    .note(:ecms_scan, "ECMS varies over the twelve scan points; the kinematic fit must use the per-point collision energy")
    .with_decay_card(decay_card_signal)
    .apply(partial_selection)

### Execute on real data, inclusive MC and the signal exclusive MC ###
root_files_full    = alg_full.execute_on(data_points + incmc_points + exMCs_signal)
root_files_partial = alg_partial.execute_on(data_points + incmc_points + exMCs_signal)