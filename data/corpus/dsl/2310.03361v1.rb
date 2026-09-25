# Paper 2310.03361v1: Cross section measurement of e+e- -> eta J/psi at sqrt(s)=3.808-4.951 GeV
# BESIII, 22.42 fb^-1 integrated luminosity, 44 energy points
# Two modes: Mode I (eta->gamma gamma) and Mode II (eta->pi0 pi+ pi-)
# J/psi reconstructed via J/psi -> l+ l- (l = e or mu)

# Representative dataset (psi(4260) peak); multi-energy scan -- see notes
scan_data = DatasetManager.real_data.find("703_4260")
scan_incMC = DatasetManager.inclusive_mc.find("703_4260")

# ---- Mode I: e+e- -> eta J/psi, J/psi -> l+l-, eta -> gamma gamma ----
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- VLL;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "ee_etajpsi_ggll_modeI_mc"
  config.related_dataset = scan_data
  config.events = 100000
  config.decay_card = decay_card_modeI
  config.cross_section = :default
end

alg_modeI = Algorithm.new("ee2etaJpsi_ModeI")
alg_modeI.set_header(["ee2etaJpsi_ModeIAlg/ee2etaJpsi_ModeI.h"])
  .set_constant({ "ECMS" => [:double, 4.260] })

sel_modeI = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot "==2"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    nlp "==1"
    nlm "==1"
  end
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_modeI.note(:lepton_separation, "Tracks with p > 1.0 GeV/c assigned as leptons; muons: E <= 0.4 GeV; electrons: E/p > 0.8")
  .note(:photon_energy_cut_modeI, "Each photon energy > 0.08 GeV after kinematic fit (radiative Bhabha suppression)")
  .note(:jpsi_mass_window, "M(l+l-) in [3.067, 3.127] GeV/c^2; sidebands [3.027,3.057] and [3.137,3.167]")
  .note(:multi_energy_scan, "This analysis runs over 44 CMS energies from 3.808 to 4.951 GeV; cross sections obtained via simultaneous unbinned ML fit of M(gamma gamma) and M(pi0 pi+ pi-) spectra")
  .note(:isr_correction, "ISR correction factor (1+delta_ISR) and vacuum polarization 1/|1-Pi|^2 applied iteratively via dressed cross section weighting method")
  .note(:helamp_generator, "Signal MC generated with HELAMP(1 0 0 0 -1 0) model")
  .note(:kinfit_chi2_cut_modeI, "chi2_4c < 40")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

alg_modeI.execute_on([scan_data, scan_incMC, exMC_modeI])

# ---- Mode II: e+e- -> eta J/psi, J/psi -> l+l-, eta -> pi0 pi+ pi- ----
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- VLL;
    Enddecay

    Decay eta
    1.000 pi0 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "ee_etajpsi_3pi2gammall_modeII_mc"
  config.related_dataset = scan_data
  config.events = 100000
  config.decay_card = decay_card_modeII
  config.cross_section = :default
end

alg_modeII = Algorithm.new("ee2etaJpsi_ModeII")
alg_modeII.set_header(["ee2etaJpsi_ModeIIAlg/ee2etaJpsi_ModeII.h"])
  .set_constant({ "ECMS" => [:double, 4.260] })

sel_modeII = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot "==4"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    identify :pion, against: [:kaon, :proton]
    nlp "==1"
    nlm "==1"
    npip "==1"
    npim "==1"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:pi0, :pip, :pim, :lp, :lm]) do
    nominal
    invariant_mass_of(:pi0, :pip, :pim).constrain_to_nominal_mass_of(:eta)
    constrain_four_momentum
    chi2_cut 200
  end

alg_modeII.note(:lepton_separation, "Tracks with p > 1.0 GeV/c assigned as leptons; muons: E <= 0.4 GeV; electrons: E/p > 0.8")
  .note(:pip_pim_separation_modeII, "Charged tracks with p <= 1.0 GeV/c assigned as pions")
  .note(:jpsi_mass_window, "M(l+l-) in [3.067, 3.127] GeV/c^2; sidebands [3.027,3.057] and [3.137,3.167]")
  .note(:multi_energy_scan, "This analysis runs over 44 CMS energies from 3.808 to 4.951 GeV")
  .note(:kinfit_chi2_cut_modeII, "chi2_5c < 80")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

alg_modeII.execute_on([scan_data, scan_incMC, exMC_modeII])