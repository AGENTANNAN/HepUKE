# arXiv: 1901.00085v3
# Observation and study of the decay J/psi -> phi eta eta'
# BESIII Collaboration
# 1.3 x 10^9 J/psi events
# Two eta' decay modes:
#   Mode I:  eta' -> gamma pi+ pi-,  eta -> gamma gamma
#   Mode II: eta' -> eta pi+ pi- (eta -> gamma gamma)

# ============================================================================
# Dataset preparation
# ============================================================================

ds_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================================
# Mode I: eta' -> gamma pi+ pi-, eta -> gamma gamma
# ============================================================================

decay_card_modeI = <<~DECAYCARD
  Decay J/psi
  1.000  phi eta eta_prime          PHSP;
  Enddecay

  Decay phi
  1.000  K+ K-                      VSS;
  Enddecay

  Decay eta
  1.000  gamma gamma                PHSP;
  Enddecay

  Decay eta_prime
  1.000  gamma pi+ pi-              PHSP;
  Enddecay

  End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Jpsi_phi_eta_etap_modeI"
  config.related_dataset = ds_jpsi
  config.events          = 500_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

algorithm_modeI = Algorithm.new("JpsiPhiEtaEtapModeI")
algorithm_modeI.set_header(["JpsiPhiEtaEtapModeIAlg/JpsiPhiEtaEtapModeI.h"])
               .set_constant("ECMS" => [:double, 3.097])

selection_modeI = Selection.new
selection_modeI
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp "==1"
    nkm "==1"
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 40
  }

algorithm_modeI
  .note(:pi0_veto, "Events with any photon pair M(γγ) in [0.12, 0.15] GeV/c^2 rejected to suppress pi0 background")
  .note(:eta_etap_assignment_modeI, "Photon assignment to eta/eta': minimize delta^2 = [M(γ1γ2)-m_eta]^2/sigma_eta^2 + [M(γ3pi+pi-)-m_etap]^2/sigma_etap^2; mass windows: eta [0.522,0.573], phi [1.010,1.030], eta' [0.936,0.979] GeV/c^2")
  .note(:mpipi_cut_modeI, "M(pi+pi-) < 0.87 GeV/c^2 to suppress J/psi -> eta phi f0(980) background")
  .note(:branching_fraction_method, "B(J/psi->phi eta eta') measured via 2D binning in M^2(phi etap) vs M^2(phi eta); efficiencies corrected area-by-area")
  .note(:x_structure, "Evidence for X structure in phi etap mass spectrum at 2.0-2.1 GeV/c^2; simultaneous fit to mode I+II; J^P=1^- significance 4.4 sigma, mass (2002.1+/-27.5+/-21.4) MeV/c^2, width (129+/-17+/-9) MeV")
  .note(:simultaneous_fit, "phi etap mass spectra fitted simultaneously across modes I and II; ROOT-level analysis")
  .with_decay_card(decay_card_modeI)
  .apply(selection_modeI)

# ============================================================================
# Mode II: eta' -> eta pi+ pi- (eta -> gamma gamma)
# ============================================================================

decay_card_modeII = <<~DECAYCARD
  Decay J/psi
  1.000  phi eta eta_prime          PHSP;
  Enddecay

  Decay phi
  1.000  K+ K-                      VSS;
  Enddecay

  Decay eta
  1.000  gamma gamma                PHSP;
  Enddecay

  Decay eta_prime
  1.000  eta pi+ pi-                PHSP;
  Enddecay

  End
DECAYCARD

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Jpsi_phi_eta_etap_modeII"
  config.related_dataset = ds_jpsi
  config.events          = 500_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

algorithm_modeII = Algorithm.new("JpsiPhiEtaEtapModeII")
algorithm_modeII.set_header(["JpsiPhiEtaEtapModeIIAlg/JpsiPhiEtaEtapModeII.h"])
                .set_constant("ECMS" => [:double, 3.097])

selection_modeII = Selection.new
selection_modeII
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp "==1"
    nkm "==1"
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 80
  }

algorithm_modeII
  .note(:pi0_veto, "Events with any photon pair M(γγ) in [0.12, 0.15] GeV/c^2 rejected to suppress pi0 background")
  .note(:eta_etap_assignment_modeII, "Photon assignment: minimize delta^2 = [M(γ1γ2)-m_eta]^2/sigma_eta^2 + [M(γ3γ4)-m_eta]^2/sigma_eta^2 for best two eta; eta whose M(pi+pi-eta) closest to m_etap assigned to eta'; mass windows: eta [0.509,0.586], phi [1.010,1.030], eta' [0.920,0.995] GeV/c^2")
  .note(:common_vertex, "All four charged tracks fitted to a common vertex before kinematic fit")
  .note(:branching_fraction_method, "B(J/psi->phi eta eta') measured via 2D binning; simultaneous fit of both modes for X structure study")
  .note(:x_structure, "Evidence for X structure in phi etap mass spectrum at 2.0-2.1 GeV/c^2; J^P=1^+ also tested (3.8 sigma)")
  .with_decay_card(decay_card_modeII)
  .apply(selection_modeII)

# ============================================================================
# Execute
# ============================================================================

root_modeI  = algorithm_modeI.execute_on([ds_jpsi, incMC_jpsi, exMC_modeI])
root_modeII = algorithm_modeII.execute_on([ds_jpsi, incMC_jpsi, exMC_modeII])