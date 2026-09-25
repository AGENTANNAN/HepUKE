# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC

# Decay card — mode I: J/psi -> gamma eta', eta' -> pi+ pi- e+ e-  (VMD model = VLL)
decay_card_etap_ee = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- e+ e- VLL;
  Enddecay

  End
DECAYCARD

# Decay card — mode II: J/psi -> gamma eta', eta' -> pi+ pi- mu+ mu-  (VMD model = VLL)
decay_card_etap_mumu = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- mu+ mu- VLL;
  Enddecay

  End
DECAYCARD

# Exclusive MC for each of the two independent final states (500k events each)
exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gametap_etap2pipiee"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_etap_ee
  config.cross_section   = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gametap_etap2pipimumu"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_etap_mumu
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# =====================================================================
# Mode I:  J/psi -> gamma eta', eta' -> pi+ pi- e+ e-
# =====================================================================
alg_name_e = "GamEtapToPiPiEE"
alg_e = Algorithm.new(alg_name_e)
alg_e.set_header(["#{alg_name_e}Alg/#{alg_name_e}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy 3.097 GeV

sel_e = Selection.new
sel_e.select_track {                 # Charged-track selection
        cos_theta 0.93               # |cos(theta)| < 0.93
        Vz        10.0               # |Vz| < 10 cm
        Vr        1.0                # Vr < 1 cm
        nChrp     "==2"              # exactly two positive tracks
        nChrn     "==2"              # exactly two negative tracks
        nNet      "==0"              # net charge zero
      }
     .select_photon {                # Photon selection
        tdc_emc_start     0          # EMC timing window 0-14
        tdc_emc_end       14
        energyThreshold_b 0.025      # barrel threshold 25 MeV
        energyThreshold_e 0.050      # endcap threshold 50 MeV
        angle_to_track    15.0       # photon-track opening angle > 15 degrees
        nGam              ">=1"      # at least one photon
      }
     .pid(method: :probability) {    # PID, probability method
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 1.0  # p>1.0 GeV -> lepton; EMC E>1.0 GeV -> e, else mu
        identify :pion, against: [:kaon]   # pi/K separation (pi+ and pi-)
        npip "==1"                   # 1 pi+
        npim "==1"                   # 1 pi-
        nlp  "==1"                   # 1 l+
        nlm  "==1"                   # 1 l-
      }
     .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {  # nominal 4C fit to gamma pi+ pi- l+ l-
        nominal                      # nominal fit — its four-momenta are stored
        constrain_four_momentum      # 4C energy-momentum constraint
        chi2_cut 200                 # loose BOSS-level chi2 cut (tight cut applied at ROOT level)
      }

# Photon-conversion veto in the e+e- mode: cannot be expressed in the DSL
# (real-photon identification from Rxy vs M_BP(e+e-) and Rxy vs Phi_ee).
alg_e.note(:background_veto,
           "photon-conversion veto applied in the e+e- mode: e+e- pairs from a "
           "converted photon are rejected using Rxy versus M_BP(e+e-) and Rxy versus "
           "Phi_ee; applied after PID and before the 4C kinematic fit")

alg_e.with_decay_card(decay_card_etap_ee).apply(sel_e)

# =====================================================================
# Mode II:  J/psi -> gamma eta', eta' -> pi+ pi- mu+ mu-
# =====================================================================
alg_name_mu = "GamEtapToPiPiMuMu"
alg_mu = Algorithm.new(alg_name_mu)
alg_mu.set_header(["#{alg_name_mu}Alg/#{alg_name_mu}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy 3.097 GeV

sel_mu = Selection.new
sel_mu.select_track {                # Charged-track selection
        cos_theta 0.93               # |cos(theta)| < 0.93
        Vz        10.0               # |Vz| < 10 cm
        Vr        1.0                # Vr < 1 cm
        nChrp     "==2"              # exactly two positive tracks
        nChrn     "==2"              # exactly two negative tracks
        nNet      "==0"              # net charge zero
      }
      .select_photon {               # Photon selection
        tdc_emc_start     0          # EMC timing window 0-14
        tdc_emc_end       14
        energyThreshold_b 0.025      # barrel threshold 25 MeV
        energyThreshold_e 0.050      # endcap threshold 50 MeV
        angle_to_track    15.0       # photon-track opening angle > 15 degrees
        nGam              ">=1"      # at least one photon
      }
      .pid(method: :probability) {   # PID, probability method
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 1.0
        identify :pion, against: [:kaon]   # pi/K separation (pi+ and pi-)
        npip "==1"                   # 1 pi+
        npim "==1"                   # 1 pi-
        nlp  "==1"                   # 1 l+
        nlm  "==1"                   # 1 l-
      }
      .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {  # nominal 4C fit to gamma pi+ pi- l+ l-
        nominal                      # nominal fit — its four-momenta are stored
        constrain_four_momentum      # 4C energy-momentum constraint
        invariant_mass_of(:lp, :lm).out_of(0.528, 0.568)  # veto |M(mu+mu-) - 0.548| < 0.02 (eta -> mu+mu-)
        chi2_cut 200                 # loose BOSS-level chi2 cut (tight chi2_4C < 25 applied at ROOT level)
      }
      # Competing-hypothesis fit: the two muons re-assigned as pions, testing
      # gamma pi+ pi- pi+ pi- (no chi2_cut, not nominal -> its chi2 is stored only).
      .assign({:lp => :pip, :lm => :pim})
      .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {
        constrain_four_momentum
      }

alg_mu.with_decay_card(decay_card_etap_mumu).apply(sel_mu)

### Execution ###
root_files_e  = alg_e.execute_on([jpsi_data, jpsi_incMC, exMC_ee])
root_files_mu = alg_mu.execute_on([jpsi_data, jpsi_incMC, exMC_mumu])