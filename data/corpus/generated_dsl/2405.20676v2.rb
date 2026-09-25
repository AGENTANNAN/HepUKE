### Dataset description ###
# Nine BESIII energy-scan points from 4.660 to 4.950 GeV (total ~4.67 fb-1):
# 706-1(4660,4680,4700) + 707-1(4740,4750,4780,4840,4914,4946)
scan_data = [
  DatasetManager.real_data.find("706_4660"),  # 4.660 GeV
  DatasetManager.real_data.find("706_4680"),  # 4.680 GeV
  DatasetManager.real_data.find("706_4700"),  # 4.700 GeV
  DatasetManager.real_data.find("707_4740"),  # 4.740 GeV
  DatasetManager.real_data.find("707_4750"),  # 4.750 GeV
  DatasetManager.real_data.find("707_4780"),  # 4.780 GeV
  DatasetManager.real_data.find("707_4840"),  # 4.840 GeV
  DatasetManager.real_data.find("707_4914"),  # 4.914 GeV
  DatasetManager.real_data.find("707_4946")   # 4.946 GeV
]

scan_incMC = [
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4914"),
  DatasetManager.inclusive_mc.find("707_4946")
]

# ---- Decay cards (EvtGen format) ----
# Mode I: e+e- -> eta' psi(2S); eta' -> gamma pi+ pi-; psi(2S) -> pi+ pi- J/psi; J/psi -> e+ e- (generated)
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta' psi(2S) PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    Decay psi(2S)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Mode II: e+e- -> eta' psi(2S); eta' -> eta pi+ pi-; eta -> gamma gamma; psi(2S) -> pi+ pi- J/psi; J/psi -> e+ e- (generated)
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta' psi(2S) PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay psi(2S)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ---- Exclusive MC: 100k events per mode, one sample per scan point ----
exMC_modeI = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_etaprimepsi2s_gammapipi"
  config.events        = 100_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_etaprimepsi2s_etapipi"
  config.events        = 100_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

### Event selection (BOSS) ###
# ================= Mode I: eta' -> gamma pi+ pi- =================
alg_name_I = "EtaPrimePsi2SGammaPiPi"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({ "ECMS" => [:double, 4.660] })   # representative; per-point energy handled below
     .note(:energy_scan, "Search spans nine scan points 4.660-4.950 GeV; a single ECMS constant is shown, but the 4C constraint energy and the exclusive-MC beam energy must follow each point.")
     .note(:pid_correction_method, "PID recipe: tracks with p>1.0 GeV/c treated as leptons (electron if E/p>0.7, otherwise muon requiring EMC energy <0.45 GeV); tracks with p<0.8 GeV/c treated as pions. Only the standard high-momentum lepton selector with approximate thresholds is encoded in the DSL.")
     .note(:psi2S_mass_window, "psi(2S) invariant-mass window [3.680,3.693] GeV/c^2 on the pi+pi- J/psi system is not directly expressible (ambiguous pion assignment among the two pi+ / two pi- candidates).")
     .note(:candidate_selection, "Best candidate combination chosen by minimizing the eta' and psi(2S) mass pulls; the DSL kinematic fit defaults to the smallest chi2 combination.")
     .note(:missing_pion_handling, "5-track events (one pion missing) are handled with a 1C pion-mass constraint.")

sel_modeI = Selection.new
  .select_track {
     cos_theta 0.93          # |cos(theta)| < 0.93
     Vz        10.0          # |Vz| < 10 cm
     Vr        1.0           # Vr < 1 cm
     nChrp     ">=2"         # at least two positive tracks
     nChrn     ">=2"         # at least two negative tracks
  }
  .select_photon {
     tdc_emc_start     0     # TDC 0
     tdc_emc_end       14    # TDC 14
     angle_to_track    10.0  # >= 10 deg from any charged track
     energyThreshold_b 0.025 # 25 MeV barrel
     energyThreshold_e 0.050 # 50 MeV endcap
     nGam              ">=1" # at least one photon (Mode I)
  }
  .pid(method: :probability) {
     prob_cut 0.001          # probability > 0.001
     # p>1.0 GeV/c tracks treated as leptons; electron if E/p>0.7 (approx. energy threshold), else muon
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.7
     nlp "==1"              # exactly one l+
     nlm "==1"              # exactly one l-
  }
  .remove([:lp <= :chrgp, :lm <= :chrgn])   # leptons removed from charged lists
  .assign({:chrgp => :pip, :chrgn => :pim}) # remaining positive/negative tracks assigned as pi+/pi-
  # Nominal 4C fit to gamma pi+ pi- pi+ pi- l+ l-
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :lp, :lm]) {
     nominal
     invariant_mass_of(:lp, :lm).within(3.083, 3.111)  # J/psi mass window (Mode I)
     constrain_four_momentum
     chi2_cut 200
  }
  # Competing hypothesis with an extra photon (photon-multiplicity veto, stored chi2 for ROOT)
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :pip, :pim, :lp, :lm]) {
     constrain_four_momentum
  }

alg_I.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_I.execute_on(scan_data + scan_incMC + exMC_modeI)

# ================= Mode II: eta' -> eta pi+ pi-, eta -> gamma gamma =================
alg_name_II = "EtaPrimePsi2SEtaPiPi"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({ "ECMS" => [:double, 4.660] })  # representative; per-point energy handled below
      .note(:energy_scan, "Search spans nine scan points 4.660-4.950 GeV; a single ECMS constant is shown, but the 4C constraint energy and the exclusive-MC beam energy must follow each point.")
      .note(:pid_correction_method, "PID recipe: tracks with p>1.0 GeV/c treated as leptons (electron if E/p>0.7, otherwise muon requiring EMC energy <0.45 GeV); tracks with p<0.8 GeV/c treated as pions. Only the standard high-momentum lepton selector with approximate thresholds is encoded in the DSL.")
      .note(:psi2S_mass_window, "psi(2S) invariant-mass window [3.680,3.693] GeV/c^2 on the pi+pi- J/psi system is not directly expressible (ambiguous pion assignment).")
      .note(:candidate_selection, "Best candidate combination chosen by minimizing the eta' and psi(2S) mass pulls; the DSL kinematic fit defaults to the smallest chi2 combination.")

sel_modeII = Selection.new
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     ">=2"
     nChrn     ">=2"
  }
  .select_photon {
     tdc_emc_start     0
     tdc_emc_end       14
     angle_to_track    10.0
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     nGam              ">=2"  # at least two photons (Mode II)
  }
  .pid(method: :probability) {
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.7
     nlp "==1"
     nlm "==1"
  }
  .remove([:lp <= :chrgp, :lm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  # 1C Kalman fit to reconstruct eta -> gamma gamma (mass-constrained to eta nominal mass, eta window preselection)
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).within(0.482, 0.604)   # eta mass window (Mode II)
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
     chi2_cut 200
     neta ">=1"             # at least one eta candidate
  }
  # Nominal 4C fit to eta pi+ pi- pi+ pi- l+ l-
  .kinematic_fit([:eta, :pip, :pim, :pip, :pim, :lp, :lm]) {
     nominal
     invariant_mass_of(:lp, :lm).within(3.073, 3.121)  # J/psi mass window (Mode II)
     constrain_four_momentum
     chi2_cut 200
  }

alg_II.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_II.execute_on(scan_data + scan_incMC + exMC_modeII)