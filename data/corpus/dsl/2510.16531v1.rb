### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")  # psi(3686) dataset at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay cards for signal: psi(3686) -> gamma chi_cJ, chi_cJ -> gamma'(X) J/psi, gamma'(X) -> e+ e-, J/psi -> l+ l-
# Three chi_cJ channels (J=0,1,2) share the same final state topology.
# Since the three chi_cJ share identical final states (gamma e+ e- l+ l-) and selection,
# a single Algorithm with the decay card covering chi_c0 is used.

# Decay card for J/psi -> e+ e- channel
decay_card_chi = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 P2GC0;
    Enddecay

    Decay chi_c0
    1.000 gamma2 J/psi PHSP;
    Enddecay

    Decay gamma2
    1.000 e+ e- PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- VLL;
    Enddecay

    End
DECAYCARD

exMC_signal_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exclusive_psi_chic_gammaX_ee_ee"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_chi
  config.cross_section = :default
end

# Decay card for J/psi -> mu+ mu- channel
decay_card_chi_mumu = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 P2GC0;
    Enddecay

    Decay chi_c0
    1.000 gamma2 J/psi PHSP;
    Enddecay

    Decay gamma2
    1.000 e+ e- PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- VLL;
    Enddecay

    End
DECAYCARD

exMC_signal_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exclusive_psi_chic_gammaX_ee_mumu"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_chi_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###
# Three chi_cJ channels (J=0,1,2) are distinguished by E_gamma and M(e+e- l+l-)
# ranges in ROOT analysis; BOSS-level selection is identical for all three.
# Final state: gamma + e+e- + e+e- (or e+e- + mu+mu-)
# 4 charged tracks (2+, 2-), at least 1 photon

# === J/psi -> e+ e- channel ===
alg_ee = Algorithm.new("PsipChicGammaX2EE")
alg_ee.set_header(["PsipChicGammaX2EEAlg/PsipChicGammaX2EE.h"])
        .set_constant({"ECMS" => [:double, 3.686]})

selection_ee = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==2"
    nlm "==2"
  }
  # Main 4C kinematic fit: constrain four-momentum
  .kinematic_fit([:gamma, :lp, :lm, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_ee
  .note(:photon_conversion_veto, "R_xy < 2 cm applied to veto photon conversion background from psi(3686)->gamma chi_cJ, chi_cJ->gamma J/psi where a photon converts to e+e- in beam pipe or MDC inner wall; applied via for_each filter or ROOT analysis")
  .note(:chic_separation, "three chi_cJ channels separated by E_gamma ranges (chi_c0: [0.23,0.28], chi_c1: [0.15,0.18], chi_c2: [0.11,0.14] GeV) and M(e+e-l+l-) ranges (chi_c0: [3.39,3.47], chi_c1: [3.48,3.55], chi_c2: [3.52,3.58] GeV); applied in ROOT analysis")
  .note(:signal_search, "search for X boson at 17 MeV/c2 and dark photon gamma' in [5,300] MeV/c2 via M(e+e-) spectrum of lower-momentum pair; unbinned extended maximum likelihood fits performed separately for each chi_cJ channel; no significant signal observed")
  .with_decay_card(decay_card_chi)
  .apply(selection_ee)

# === J/psi -> mu+ mu- channel ===
alg_mumu = Algorithm.new("PsipChicGammaX2MuMu")
alg_mumu.set_header(["PsipChicGammaX2MuMuAlg/PsipChicGammaX2MuMu.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

selection_mumu = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.4
    nlp "==2"
    nlm "==2"
  }
  .kinematic_fit([:gamma, :lp, :lm, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mumu
  .note(:photon_conversion_veto, "R_xy < 2 cm applied to veto photon conversion background; applied via for_each filter or ROOT analysis")
  .note(:chic_separation, "three chi_cJ channels separated by E_gamma and M(e+e-mu+mu-) mass ranges in ROOT analysis")
  .note(:signal_search, "search for X boson at 17 MeV/c2 and dark photon gamma' in [5,300] MeV/c2; unblinded fits gave eps_c < 1.2e-2 at 90% CL and dark photon mixing epsilon < (2.5-17.5)e-3")
  .with_decay_card(decay_card_chi_mumu)
  .apply(selection_mumu)

# Execute the algorithms
root_files_ee = alg_ee.execute_on([psip_data, psip_incMC, exMC_signal_ee])
root_files_mumu = alg_mumu.execute_on([psip_data, psip_incMC, exMC_signal_mumu])