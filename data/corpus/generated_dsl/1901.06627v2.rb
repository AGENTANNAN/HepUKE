### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) real data sample
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card: psi(2S) -> gamma chi_cJ, chi_cJ -> mu+ mu- J/psi, J/psi -> e+ e-
decay_card_ee = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.000 mu+ mu- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card: psi(2S) -> gamma chi_cJ, chi_cJ -> mu+ mu- J/psi, J/psi -> mu+ mu-
decay_card_mumu = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.000 mu+ mu- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC sample for each J/psi decay mode
exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachi_mumu_ee"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_ee
  config.cross_section   = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachi_mumu_mumu"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Common charged-track and photon selection shared by both J/psi decay modes
common_selection = Selection.new
    .select_track {
        cos_theta 0.93       # |cos(theta)| < 0.93
        Vz        10.0       # |Vz| < 10 cm
        Vr        1.0        # Vr < 1 cm
        nChrp     "==2"      # exactly two positive tracks
        nChrn     "==2"      # exactly two negative tracks
        nNet      "==0"      # net charge zero
    }
    .select_photon {
        tdc_emc_start     0      # EMC timing window 0-14
        tdc_emc_end       14
        angle_to_track    20.0   # angle to any charged track > 20 degrees
        energyThreshold_b 0.025  # 25 MeV barrel threshold
        energyThreshold_e 0.050  # 50 MeV endcap threshold
        nGam              ">=1"  # at least one photon
    }

# ---- Mode I: J/psi -> e+ e- ----
selection_ee = common_selection.dup
    .pid(method: :probability) {
        # high-momentum tracks are treated as leptons; a lepton is an electron if
        # its EMC energy exceeds 1.0 GeV, otherwise a muon
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 1.0
        # the mu+ mu- pair from chi_cJ is identified against pi/K/p
        identify :muon, against: [:pion, :kaon, :proton]
        nlp  "==1"      # exactly one l+
        nlm  "==1"      # exactly one l-
        nmup ">=1"      # at least one mu+
        nmum ">=1"      # at least one mu-
    }
    .kinematic_fit([:gamma, :mup, :mum, :lp, :lm]) {
        nominal
        vertex_fit([1, 2, 3, 4])                           # vertex fit on the four charged tracks
        constrain_four_momentum                            # 4C constraint
        invariant_mass_of(:lp, :lm).between(3.085, 3.110)  # J/psi lepton-pair mass window (GeV/c^2)
        invariant_mass_of(:gamma, :mup, :mum).out_of(0.535, 0.560)  # suppress psi(2S) -> eta J/psi
        chi2_cut 40
    }

# ---- Mode II: J/psi -> mu+ mu- ----
selection_mumu = common_selection.dup
    .pid(method: :probability) {
        # high-momentum tracks are treated as leptons; a lepton is an electron if
        # its EMC energy exceeds 0.3 GeV, otherwise a muon
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.3
        identify :muon, against: [:pion, :kaon, :proton]
        nmup ">=2"      # at least two mu+
        nmum ">=2"      # at least two mu-
    }
    .kinematic_fit([:gamma, :mup, :mum, :lp, :lm]) {
        nominal
        vertex_fit([1, 2, 3, 4])                           # vertex fit on the four charged tracks
        constrain_four_momentum                            # 4C constraint
        invariant_mass_of(:lp, :lm).between(3.085, 3.110)  # J/psi lepton-pair mass window (GeV/c^2)
        invariant_mass_of(:gamma, :mup, :mum).out_of(0.535, 0.560)  # suppress psi(2S) -> eta J/psi
        chi2_cut 40
    }

### Algorithms — one per independent J/psi decay mode (different PID and fit hypotheses) ###
alg_ee = Algorithm.new("GamChiCJpsiEE")
alg_ee.set_header(["GamChiCJpsiEEAlg/GamChiCJpsiEE.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
      .note(:helix_correction, "helix-parameter correction of the charged tracks, obtained from the psi(3686) -> pi+ pi- J/psi control sample, is applied before the 4C kinematic fit; the corrected helix parameters are used as inputs to the kinematic fit")
      .note(:background_veto, "R_xy > 8.5 cm requirement on the photon candidate is used to veto photon conversions; the conversion-point radius is not available through the standard selection methods and is applied as an external step before the kinematic fit")

alg_mumu = Algorithm.new("GamChiCJpsiMuMu")
alg_mumu.set_header(["GamChiCJpsiMuMuAlg/GamChiCJpsiMuMu.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .note(:helix_correction, "helix-parameter correction of the charged tracks, obtained from the psi(3686) -> pi+ pi- J/psi control sample, is applied before the 4C kinematic fit; the corrected helix parameters are used as inputs to the kinematic fit")
        .note(:background_veto, "R_xy > 8.5 cm requirement on the photon candidate is used to veto photon conversions; the conversion-point radius is not available through the standard selection methods and is applied as an external step before the kinematic fit")

# Attach the corresponding decay card and render the selection into BOSS code.
# When several photon candidates are present, the kinematic fit automatically
# keeps the combination giving the smallest fit chi2.
alg_ee.with_decay_card(decay_card_ee).apply(selection_ee)
alg_mumu.with_decay_card(decay_card_mumu).apply(selection_mumu)

# Execute on real data, inclusive MC and the mode-specific exclusive MC
root_files_ee   = alg_ee.execute_on([psip_data, psip_incMC, exMC_ee])
root_files_mumu = alg_mumu.execute_on([psip_data, psip_incMC, exMC_mumu])