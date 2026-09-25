# =====================================================================
#  e+e- -> eta' J/psi  @ 14 energy points (sqrt(s) = 4.178 - 4.600 GeV)
#  BOSS part: dataset preparation + event selection (up to the final fit)
# =====================================================================

### Dataset description ###
# 14 real-data points (BOSS 7.0.3); sample name = "<BOSS>_<CMS energy in MeV>".
energy_points = [
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4600")
]

# Corresponding inclusive MC sample at each energy point.
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4237"),
  DatasetManager.inclusive_mc.find("703_4246"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4270"),
  DatasetManager.inclusive_mc.find("703_4280"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4600")
]

### Decay cards (EvtGen format) ###
# Mode I: psi(4260) -> eta' J/psi, eta' -> pi+ pi- eta, eta -> gamma gamma, J/psi -> e+ e-
decay_card_modeI = <<~DECAYCARD
  Decay psi(4260)
  1.000 eta' J/psi PHSP;
  Enddecay

  Decay eta'
  1.000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Mode II: psi(4260) -> eta' J/psi, eta' -> gamma pi+ pi-, J/psi -> e+ e-
decay_card_modeII = <<~DECAYCARD
  Decay psi(4260)
  1.000 eta' J/psi PHSP;
  Enddecay

  Decay eta'
  1.000 gamma pi+ pi- PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

### Exclusive MC: 100,000 events at each energy point, for each decay mode ###
exMC_modeI = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_etapJpsi_modeI_pipieta"
  config.events        = 100_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_etapJpsi_modeII_gampipi"
  config.events        = 100_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ------------------------------ Mode I ------------------------------
# eta' -> pi+ pi- eta, eta -> gamma gamma, J/psi -> l+ l-
alg_name_modeI = "EtapJpsiModeI"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:helix_correction, "charged-track helix parameters are corrected before the kinematic fits")
         .note(:background_veto, "photon conversion suppressed by requiring cos(theta(pi+,pi-)) < 0.95")
         .note(:pid_correction_method, "high-momentum lepton PID: p>1.0 GeV treated as a lepton, EMC energy >1.0 GeV as an electron and <0.4 GeV as a muon (the 0.4 GeV muon threshold is not separately tunable through identify_high_momentum_leptons)")

sel_modeI = Selection.new
sel_modeI.select_track {
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     "==2"
            nChrn     "==2"
            nNet      "==0"
          }
          .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    20.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=2"    # Mode I needs at least two photons
          }
          .pid(method: :probability) {
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 1.0
            identify :pion, against: [:kaon, :proton]
            nlp   "==1"
            nlm   "==1"
            npip  "==1"
            npim  "==1"
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
            chi2_cut 25
            neta     ">=1"
          }
          .kinematic_fit([:lp, :lm, :pip, :pim, :eta]) {
            nominal                  # 5C = 4C (energy-momentum) here + the 1C eta mass constraint from the Kalman fit
            constrain_four_momentum
            invariant_mass_of(:lp, :lm).within(3.07, 3.13)              # J/psi window (sidebands handled in ROOT)
            invariant_mass_of(:eta, :lp, :lm).out_of(3.67, 3.70)        # veto psi(2S) -> eta J/psi
            chi2_cut 50
          }   # combination with the smallest chi2 is chosen by default

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
root_files_modeI = alg_modeI.execute_on(energy_points + incMC_points + exMC_modeI)

# ------------------------------ Mode II -----------------------------
# eta' -> gamma pi+ pi-, J/psi -> l+ l-
alg_name_modeII = "EtapJpsiModeII"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:helix_correction, "charged-track helix parameters are corrected before the kinematic fit")
          .note(:background_veto, "photon conversion suppressed by requiring cos(theta(pi+,pi-)) < 0.95")
          .note(:pid_correction_method, "high-momentum lepton PID: p>1.0 GeV treated as a lepton, EMC energy >1.0 GeV as an electron and <0.4 GeV as a muon (the 0.4 GeV muon threshold is not separately tunable through identify_high_momentum_leptons)")

sel_modeII = Selection.new
sel_modeII.select_track {
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     "==2"
            nChrn     "==2"
            nNet      "==0"
          }
          .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    20.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=1"    # Mode II needs at least one photon
          }
          .pid(method: :probability) {
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 1.0
            identify :pion, against: [:kaon, :proton]
            nlp   "==1"
            nlm   "==1"
            npip  "==1"
            npim  "==1"
          }
          .kinematic_fit([:gamma, :lp, :lm, :pip, :pim]) {
            nominal
            constrain_four_momentum
            invariant_mass_of(:lp, :lm).within(3.07, 3.13)              # J/psi window (sidebands handled in ROOT)
            invariant_mass_of(:pip, :pim, :lp, :lm).out_of(3.66, 3.71)  # veto psi(2S) -> pi+ pi- J/psi
            chi2_cut 40
          }   # combination with the smallest chi2 is chosen by default

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
root_files_modeII = alg_modeII.execute_on(energy_points + incMC_points + exMC_modeII)