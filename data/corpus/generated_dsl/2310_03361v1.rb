# Core DSL classes and dependencies are loaded automatically at execution

### Dataset preparation ###
# 44 e+e- energy points in the range 3.808 - 4.951 GeV (BOSS_energy sample names)
sample_names = %w[
  703_3810 703_3900 703_4009 703_4090 703_4180 703_4190 703_4200 703_4210
  703_4220 703_4230 703_4237 703_4245 703_4246 703_4260 703_4270 703_4280
  703_4310 703_4360 703_4390 703_4420 703_4470 703_4530 703_4575 703_4600
  705_4130 705_4160 705_4290 705_4315 705_4340 705_4380 705_4400 705_4440
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]

# Corresponding real data and inclusive MC at all energy points
data_points  = sample_names.map { |s| DatasetManager.real_data.find(s) }
incMC_points = sample_names.map { |s| DatasetManager.inclusive_mc.find(s) }

# Shared decay card: e+e- -> eta J/psi, with eta -> gamma gamma (Mode I) /
# eta -> pi0 pi+ pi- (Mode II), and J/psi -> e+e- / mu+mu-
decay_card_etaJpsi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta J/psi PHSP;
  Enddecay

  Decay eta
  0.3941 gamma gamma PHSP;
  0.2274 pi0 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  0.0594 e+ e- PHOTOS VLL;
  0.0593 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive MC for Mode I, generated at every scan point with the shared card
exMC_modeI = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name  = "sig_etaJpsi_modeI_ggll"
  config.events       = 100_000
  config.decay_card   = decay_card_etaJpsi
  config.cross_section = :default
end

# 100k-event exclusive MC for Mode II, same e+e- -> eta J/psi decay card
exMC_modeII = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name  = "sig_etaJpsi_modeII_pi0pipimll"
  config.events       = 100_000
  config.decay_card   = decay_card_etaJpsi
  config.cross_section = :default
end

### Event selection (BOSS) — Mode I: e+e- -> eta J/psi, eta -> gamma gamma, J/psi -> l+l- ###
alg_name_I = "EtaJpsiModeI"
alg_modeI = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})   # nominal scan energy (per-point handled at runtime)
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:lepton_id_recipe,
               "high-momentum tracks (p>1.0 GeV/c) treated as leptons; electron-like if E/pc>0.8, " \
               "muon-like if E<=0.4 GeV. identify_high_momentum_leptons approximates the e/mu split " \
               "with an ERaw threshold; at most one same-flavour opposite-charge pair (e+e- or mu+mu-) " \
               "is accepted as the J/psi candidate.")

sel_modeI = Selection.new
sel_modeI.select_track {                 # two charged tracks (l+ l-)
           cos_theta 0.93                # |cos(theta)| < 0.93
           Vz        10.0                # |Vz| < 10 cm
           Vr        1.0                 # Vr < 1 cm
           nChrp     "==1"               # exactly 2 tracks in total
           nChrn     "==1"
           nNet      "==0"               # zero net charge
         }
         .select_photon {                # at least two photons
           tdc_emc_start     0           # EMC timing window
           tdc_emc_end       14
           energyThreshold_b 0.025       # barrel threshold 25 MeV
           energyThreshold_e 0.050       # endcap threshold 50 MeV
           angle_to_track    20.0        # opening angle to any track > 20 deg
           nGam              ">=2"
         }
         .pid(method: :probability) {    # probability PID
           prob_cut 0.001
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
           nlp "==1"                     # one positive lepton
           nlm "==1"                     # one negative lepton (opposite charge)
         }
         .kinematic_fit([:gamma, :gamma, :lp, :lm]) {   # 4C fit to gamma gamma l+ l-
           nominal
           constrain_four_momentum
           chi2_cut 40
         }

alg_modeI.with_decay_card(decay_card_etaJpsi).apply(sel_modeI)
alg_modeI.execute_on(data_points + incMC_points + exMC_modeI)

### Event selection (BOSS) — Mode II: e+e- -> eta J/psi, eta -> pi0 pi+ pi-, J/psi -> l+l- ###
alg_name_II = "EtaJpsiModeII"
alg_modeII = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:lepton_id_recipe,
                "high-momentum tracks (p>1.0 GeV/c) treated as leptons; electron-like if E/pc>0.8, " \
                "muon-like if E<=0.4 GeV; only same-flavour opposite-charge lepton pairs are accepted " \
                "as the J/psi candidate.")

sel_modeII = Selection.new
sel_modeII.select_track {                # four charged tracks (pi+ pi- l+ l-)
            cos_theta 0.93               # |cos(theta)| < 0.93
            Vz        10.0               # |Vz| < 10 cm
            Vr        1.0                # Vr < 1 cm
            nChrp     "==2"              # exactly 4 tracks in total
            nChrn     "==2"
            nNet      "==0"              # zero net charge
          }
          .select_photon {               # at least two photons (from pi0)
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track    20.0       # opening angle to any track > 20 deg
            nGam              ">=2"
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 0.6
            identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K and p
            nlp  "==1"
            nlm  "==1"
            npip "==1"
            npim "==1"
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {     # 1C pi0 reconstruction from two photons
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
          }
          .kinematic_fit([:gamma, :gamma, :pip, :pim, :lp, :lm]) {   # 5C fit: 4-momentum + pi0 mass
            nominal
            constrain_four_momentum
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 80
          }

alg_modeII.with_decay_card(decay_card_etaJpsi).apply(sel_modeII)
alg_modeII.execute_on(data_points + incMC_points + exMC_modeII)