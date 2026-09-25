# DSL for arxiv:1812.05800v1
# psi(3686) -> p pbar eta' and J/psi -> p pbar eta'
# Two eta' decay modes: eta' -> gamma pi+ pi- and eta' -> eta pi+ pi- (eta -> gamma gamma)
# 4 separate Algorithm objects: 2 modes x 2 datasets

# Common decay card for psi(3686) -> p pbar eta'
decay_card_psip = <<~DECAYCARD
  Decay psi(100443)
  1 p+ anti-p- eta' PHSP;
  Enddecay
  Decay eta'
  1 gamma pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

decay_card_psip_eta = <<~DECAYCARD
  Decay psi(100443)
  1 p+ anti-p- eta' PHSP;
  Enddecay
  Decay eta'
  1 eta pi+ pi- PHSP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Common decay card for J/psi -> p pbar eta'
decay_card_jpsi = <<~DECAYCARD
  Decay J/psi
  1 p+ anti-p- eta' PHSP;
  Enddecay
  Decay eta'
  1 gamma pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

decay_card_jpsi_eta = <<~DECAYCARD
  Decay J/psi
  1 p+ anti-p- eta' PHSP;
  Enddecay
  Decay eta'
  1 eta pi+ pi- PHSP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Datasets
data_psip = DatasetManager.real_data.find("709_3686")
incMC_psip = DatasetManager.inclusive_mc.find("709_3686")
data_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# =====================
# Mode 1: eta' -> gamma pi+ pi-
# =====================

# -- psi(3686) --
sig_psip_gamma_pipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_ppbar_etap_gamma_pipi"
  config.related_dataset = data_psip
  config.events          = 100_000
  config.decay_card      = decay_card_psip
  config.cross_section   = :default
end

alg_psip_gamma_pipi = Algorithm.new("PsipEtaPGamPiPi")
alg_psip_gamma_pipi.set_header(["PsipEtaPGamPiPiAlg/PsipEtaPGamPiPi.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })

sel_psip_gamma_pipi = Selection.new
  .select_track do
    nChrp "==0"
    nChrn "==0"
    nTot "==4"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .assign(:prp, from: :positive)
  .assign(:prm, from: :negative)
  .select_photon do
    nPhoton ">=1"
    photon_energy 0.025
    min_angle_to_charged 10.0
  end
  .pid do
    prob_cut 0.001
    identify :pip, "from_kaon_and_proton"
    identify :pim, "from_kaon_and_proton"
    identify :prp, "from_pion_and_kaon"
    identify :prm, "from_pion_and_kaon"
  end
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_psip_gamma_pipi.with_decay_card(decay_card_psip).apply(sel_psip_gamma_pipi)
alg_psip_gamma_pipi.note(:invariant_mass_cut, "920 < M(gamma pi+ pi-) < 1000 MeV/c^2 (eta' signal region), sideband 932.4-942.0 and 974.0-983.6 MeV/c^2")
alg_psip_gamma_pipi.note(:chi_cj_veto, "Veto psi(3686) -> gamma chi_cJ, chi_cJ -> p pbar pi+ pi- by recoil mass windows")
alg_psip_gamma_pipi.note(:jpsi_veto, "Veto psi(3686) -> pi+ pi- J/psi, J/psi -> gamma p pbar by recoil mass windows")
alg_psip_gamma_pipi.note(:eta_sideband_fit, "Signal yield from simultaneous unbinned ML fit to gamma pi+ pi- and eta pi+ pi- invariant mass spectra")
alg_psip_gamma_pipi.note(:photon_multiple, "When more photons than required, loop over combinations and keep minimal chi2 from kinematic fit")
alg_psip_gamma_pipi.note(:eta_prime_window, "eta' signal region: 948.2 < M < 967.4 MeV/c^2")
alg_psip_gamma_pipi.execute_on([data_psip, incMC_psip, sig_psip_gamma_pipi])

# -- J/psi --
sig_jpsi_gamma_pipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_ppbar_etap_gamma_pipi"
  config.related_dataset = data_jpsi
  config.events          = 100_000
  config.decay_card      = decay_card_jpsi
  config.cross_section   = :default
end

alg_jpsi_gamma_pipi = Algorithm.new("JpsiEtaPGamPiPi")
alg_jpsi_gamma_pipi.set_header(["JpsiEtaPGamPiPiAlg/JpsiEtaPGamPiPi.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_jpsi_gamma_pipi = Selection.new
  .select_track do
    nChrp "==0"
    nChrn "==0"
    nTot "==4"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .assign(:prp, from: :positive)
  .assign(:prm, from: :negative)
  .select_photon do
    nPhoton ">=1"
    photon_energy 0.025
    min_angle_to_charged 10.0
  end
  .pid do
    prob_cut 0.001
    identify :pip, "from_kaon_and_proton"
    identify :pim, "from_kaon_and_proton"
    identify :prp, "from_pion_and_kaon"
    identify :prm, "from_pion_and_kaon"
  end
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_jpsi_gamma_pipi.with_decay_card(decay_card_jpsi).apply(sel_jpsi_gamma_pipi)
alg_jpsi_gamma_pipi.note(:invariant_mass_cut, "920 < M(gamma pi+ pi-) < 1000 MeV/c^2; eta' signal region: 948.2 < M < 967.4 MeV/c^2")
alg_jpsi_gamma_pipi.note(:eta_sideband_fit, "Signal yield from simultaneous unbinned ML fit to gamma pi+ pi- and eta pi+ pi- invariant mass spectra")
alg_jpsi_gamma_pipi.note(:photon_multiple, "When more photons than required, loop over combinations and keep minimal chi2")
alg_jpsi_gamma_pipi.execute_on([data_jpsi, incMC_jpsi, sig_jpsi_gamma_pipi])

# =====================
# Mode 2: eta' -> eta pi+ pi- (eta -> gamma gamma)
# =====================

# -- psi(3686) --
sig_psip_eta_pipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_ppbar_etap_eta_pipi"
  config.related_dataset = data_psip
  config.events          = 100_000
  config.decay_card      = decay_card_psip_eta
  config.cross_section   = :default
end

alg_psip_eta_pipi = Algorithm.new("PsipEtaPEtaPiPi")
alg_psip_eta_pipi.set_header(["PsipEtaPEtaPiPiAlg/PsipEtaPEtaPiPi.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })

sel_psip_eta_pipi = Selection.new
  .select_track do
    nChrp "==0"
    nChrn "==0"
    nTot "==4"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .assign(:prp, from: :positive)
  .assign(:prm, from: :negative)
  .select_photon do
    nPhoton ">=2"
    photon_energy 0.025
    min_angle_to_charged 10.0
  end
  .pid do
    prob_cut 0.001
    identify :pip, "from_kaon_and_proton"
    identify :pim, "from_kaon_and_proton"
    identify :prp, "from_pion_and_kaon"
    identify :prm, "from_pion_and_kaon"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    constrain_to_nominal_mass_of :eta
  end
  .kinematic_fit([:prp, :prm, :pip, :pim, :eta]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_psip_eta_pipi.with_decay_card(decay_card_psip_eta).apply(sel_psip_eta_pipi)
alg_psip_eta_pipi.note(:invariant_mass_cut, "920 < M(eta pi+ pi-) < 1000 MeV/c^2; eta' signal region: 948.2 < M < 967.4 MeV/c^2")
alg_psip_eta_pipi.note(:jpsi_veto_etapip, "Veto psi(3686) -> eta J/psi, eta->gamma gamma, J/psi -> p pbar pi+ pi- by recoil mass windows")
alg_psip_eta_pipi.note(:jpsi_veto_pipi, "Veto psi(3686) -> pi+ pi- J/psi, J/psi -> eta p pbar by recoil mass windows")
alg_psip_eta_pipi.note(:eta_sideband_fit, "Signal yield from simultaneous unbinned ML fit")
alg_psip_eta_pipi.note(:photon_multiple, "When more photons than required, loop over combinations and keep minimal chi2")
alg_psip_eta_pipi.execute_on([data_psip, incMC_psip, sig_psip_eta_pipi])

# -- J/psi --
sig_jpsi_eta_pipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_ppbar_etap_eta_pipi"
  config.related_dataset = data_jpsi
  config.events          = 100_000
  config.decay_card      = decay_card_jpsi_eta
  config.cross_section   = :default
end

alg_jpsi_eta_pipi = Algorithm.new("JpsiEtaPEtaPiPi")
alg_jpsi_eta_pipi.set_header(["JpsiEtaPEtaPiPiAlg/JpsiEtaPEtaPiPi.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_jpsi_eta_pipi = Selection.new
  .select_track do
    nChrp "==0"
    nChrn "==0"
    nTot "==4"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .assign(:prp, from: :positive)
  .assign(:prm, from: :negative)
  .select_photon do
    nPhoton ">=2"
    photon_energy 0.025
    min_angle_to_charged 10.0
  end
  .pid do
    prob_cut 0.001
    identify :pip, "from_kaon_and_proton"
    identify :pim, "from_kaon_and_proton"
    identify :prp, "from_pion_and_kaon"
    identify :prm, "from_pion_and_kaon"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    constrain_to_nominal_mass_of :eta
  end
  .kinematic_fit([:prp, :prm, :pip, :pim, :eta]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_jpsi_eta_pipi.with_decay_card(decay_card_jpsi_eta).apply(sel_jpsi_eta_pipi)
alg_jpsi_eta_pipi.note(:invariant_mass_cut, "920 < M(eta pi+ pi-) < 1000 MeV/c^2; eta' signal region: 948.2 < M < 967.4 MeV/c^2")
alg_jpsi_eta_pipi.note(:eta_sideband_fit, "Signal yield from simultaneous unbinned ML fit")
alg_jpsi_eta_pipi.note(:photon_multiple, "When more photons than required, loop over combinations and keep minimal chi2")
alg_jpsi_eta_pipi.execute_on([data_jpsi, incMC_jpsi, sig_jpsi_eta_pipi])