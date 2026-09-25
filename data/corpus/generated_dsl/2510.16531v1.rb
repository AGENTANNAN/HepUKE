# Core DSL classes and dependencies are loaded automatically at execution

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) real data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # matching inclusive MC sample

# Decay card (written for chi_c0; chi_c1 and chi_c2 share the identical topology)
# --- e+ e- channel:  J/psi -> e+ e- ---
decay_card_ee = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0           P2GC1;
    Enddecay

    Decay chi_c0
    1.0000 X J/psi                PHSP;
    Enddecay

    Decay X
    1.0000 e+ e-                  PHOTOS VLL;
    Enddecay

    Decay J/psi
    1.0000 e+ e-                  PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# --- mu+ mu- channel:  J/psi -> mu+ mu- ---
decay_card_mumu = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0           P2GC1;
    Enddecay

    Decay chi_c0
    1.0000 X J/psi                PHSP;
    Enddecay

    Decay X
    1.0000 e+ e-                  PHOTOS VLL;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu-                PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for each J/psi lepton mode
exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammachic0_ee"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_ee
  config.cross_section   = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammachic0_mumu"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Common charged-track and photon selection shared by both lepton channels
event_selection_common = Selection.new
    .select_track {
        cos_theta 0.93        # |cos(theta)| < 0.93
        Vz        10.0        # |Vz| < 10 cm
        Vr        1.0         # Vr < 1 cm
        nChrp     "==2"       # exactly two positive tracks
        nChrn     "==2"       # exactly two negative tracks
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0    # angle to nearest charged track > 10 deg
        energyThreshold_b 0.025   # 25 MeV in the barrel
        energyThreshold_e 0.050   # 50 MeV in the endcap
        nGam   ">=1"              # at least one photon
    }

# ---------------- e+ e- channel (electron threshold 0.6 GeV) ----------------
alg_name_ee = "GamChiC0EE"
alg_ee = Algorithm.new(alg_name_ee)
alg_ee.set_header(["#{alg_name_ee}Alg/#{alg_name_ee}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

selection_ee = event_selection_common.dup
    .pid(method: :probability) {
        # high-momentum tracks (p > 1.0 GeV) treated as leptons; electron if EMC eraw > 0.6 GeV, else muon
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                        treat_as_electron_if_energy_above: 0.6
        nlp "==2"   # two l+ (here e+)
        nlm "==2"   # two l- (here e-)
    }
    .kinematic_fit([:gamma, :lp, :lm, :lp, :lm]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

alg_ee.note(:background_veto, "photon conversions vetoed by requiring R_xy < 2 cm on the reconstructed lepton tracks")
alg_ee.with_decay_card(decay_card_ee).apply(selection_ee)
root_files_ee = alg_ee.execute_on([psip_data, psip_incMC, exMC_ee])

# ---------------- mu+ mu- channel (electron threshold 0.4 GeV) ----------------
alg_name_mumu = "GamChiC0MuMu"
alg_mumu = Algorithm.new(alg_name_mumu)
alg_mumu.set_header(["#{alg_name_mumu}Alg/#{alg_name_mumu}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})

selection_mumu = event_selection_common.dup
    .pid(method: :probability) {
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                        treat_as_electron_if_energy_above: 0.4
        nlp "==2"   # two l+ (here e+ and mu+)
        nlm "==2"   # two l- (here e- and mu-)
    }
    .kinematic_fit([:gamma, :lp, :lm, :lp, :lm]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

alg_mumu.note(:background_veto, "photon conversions vetoed by requiring R_xy < 2 cm on the reconstructed lepton tracks")
alg_mumu.with_decay_card(decay_card_mumu).apply(selection_mumu)
root_files_mumu = alg_mumu.execute_on([psip_data, psip_incMC, exMC_mumu])