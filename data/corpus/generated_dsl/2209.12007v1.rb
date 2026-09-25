# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
# e+e- -> pi+ pi- J/psi is studied at four c.m. energies:
#   3807.7 MeV -> sample 703_3810
#   3867.4 MeV -> sample 703_3872
#   3871.3 MeV -> sample 703_3871
#   3896.2 MeV -> sample 703_3900
data_points = [
  DatasetManager.real_data.find("703_3810"),   # 3807.7 MeV
  DatasetManager.real_data.find("703_3872"),   # 3867.4 MeV
  DatasetManager.real_data.find("703_3871"),   # 3871.3 MeV
  DatasetManager.real_data.find("703_3900")    # 3896.2 MeV
]
incmc_points = [
  DatasetManager.inclusive_mc.find("703_3810"),
  DatasetManager.inclusive_mc.find("703_3872"),
  DatasetManager.inclusive_mc.find("703_3871"),
  DatasetManager.inclusive_mc.find("703_3900")
]

# Decay card for e+e- -> pi+ pi- J/psi, J/psi -> e+ e-
decay_card_ee = <<~DECAYCARD
  Decay psi(4260)
  1  pi+  pi-  J/psi    PHSP;
  Enddecay

  Decay J/psi
  1  e+  e-    PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Decay card for e+e- -> pi+ pi- J/psi, J/psi -> mu+ mu-
decay_card_mumu = <<~DECAYCARD
  Decay psi(4260)
  1  pi+  pi-  J/psi    PHSP;
  Enddecay

  Decay J/psi
  1  mu+  mu-    PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# 500k-event exclusive MC per leptonic mode, generated at every energy point
exMC_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_pipiJpsi_ee"
  config.events        = 500_000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_pipiJpsi_mumu"
  config.events        = 500_000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ---------- Channel I : J/psi -> e+ e- ----------
alg_name_ee = "PipiJpsiEE"
alg_ee = Algorithm.new(alg_name_ee)
alg_ee.set_header(["#{alg_name_ee}Alg/#{alg_name_ee}.h"])
      .set_constant({ "ECMS" => [:double, 3.896] })
      .set_alias({ "std::vector<double>" => "Vdouble" })

sel_ee = Selection.new
  .select_track {                       # four charged tracks: |cos(theta)|<0.93, |Vz|<10 cm, Vr<1 cm
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"                     # net charge zero
  }
  .pid(method: :probability) {          # probability PID, prob > 0.001
    prob_cut 0.001
    identify :electron, against: [:pion]            # e+/e- identified against pions
    identify :pion, against: [:electron, :kaon]     # pi+/pi- identified against electrons and kaons
    nep  "==1"                          # one e+
    nem  "==1"                          # one e-
    npip "==1"                          # one pi+
    npim "==1"                          # one pi-
  }
  .remove(:pip) { condition "three_momentum_of(:pip) > 0.6" }   # pion momentum < 0.6 GeV/c
  .remove(:pim) { condition "three_momentum_of(:pim) > 0.6" }
  .remove(:ep)  { condition "three_momentum_of(:ep) < 1.0" }    # lepton momentum > 1.0 GeV/c
  .remove(:em)  { condition "three_momentum_of(:em) < 1.0" }
  .kinematic_fit([:pip, :pim, :ep, :em]) {   # 4C kinematic fit, nominal hypothesis
    nominal
    constrain_four_momentum
    chi2_cut 200                        # loose cut in BOSS; tight chi2 < 60 applied in ROOT
  }

alg_ee.note(:electron_ep_cut,
            "electron candidates required to satisfy E/p > 1.1 (EMC shower energy over MDC momentum); " \
            "no dedicated DSL primitive available")
     .note(:background_veto,
            "photon-conversion vetoes: reject events with cos(theta)(pi+pi-) > 0.95 or " \
            "cos(theta)(pi+- e-+) > 0.98; angular-correlation cuts not expressible in the current DSL")
     .with_decay_card(decay_card_ee)
     .apply(sel_ee)

# ---------- Channel II : J/psi -> mu+ mu- ----------
alg_name_mumu = "PipiJpsiMuMu"
alg_mumu = Algorithm.new(alg_name_mumu)
alg_mumu.set_header(["#{alg_name_mumu}Alg/#{alg_name_mumu}.h"])
        .set_constant({ "ECMS" => [:double, 3.896] })
        .set_alias({ "std::vector<double>" => "Vdouble" })

sel_mumu = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :muon, against: [:pion]     # mu+/mu- identified against pions
    identify :pion, against: [:kaon]     # pi+/pi- identified against kaons
    nmup "==1"                           # one mu+
    nmum "==1"                           # one mu-
    npip "==1"                           # one pi+
    npim "==1"                           # one pi-
  }
  .remove(:pip) { condition "three_momentum_of(:pip) > 0.6" }   # pion momentum < 0.6 GeV/c
  .remove(:pim) { condition "three_momentum_of(:pim) > 0.6" }
  .remove(:mup) { condition "three_momentum_of(:mup) < 1.0" }   # lepton momentum > 1.0 GeV/c
  .remove(:mum) { condition "three_momentum_of(:mum) < 1.0" }
  .kinematic_fit([:pip, :pim, :mup, :mum]) {  # 4C kinematic fit, nominal hypothesis
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mumu.note(:muon_ep_cut,
              "muon candidates required to satisfy E/p < 0.35 (EMC shower energy over MDC momentum); " \
              "no dedicated DSL primitive available")
       .note(:background_veto,
              "photon-conversion vetoes: reject events with cos(theta)(pi+pi-) > 0.95 or " \
              "cos(theta)(pi+- mu-+) > 0.98; angular-correlation cuts not expressible in the current DSL")
       .with_decay_card(decay_card_mumu)
       .apply(sel_mumu)

# Execute both channels on real data, inclusive MC and the per-energy exclusive MC
root_files_ee   = alg_ee.execute_on(data_points + incmc_points + exMC_ee)
root_files_mumu = alg_mumu.execute_on(data_points + incmc_points + exMC_mumu)