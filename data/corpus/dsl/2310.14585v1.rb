# psi(3686) -> Sigma+ anti-Sigma- omega and Sigma+ anti-Sigma- phi
# BESIII, arXiv:2310.14585v1
# First observations of these decay modes
# Two independent decay modes: omega-mode and phi-mode

### Dataset description ###

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ============================================================
# Mode I: psi(3686) -> Sigma+ anti-Sigma- omega
# Sigma+ -> p pi0, anti-Sigma- -> anti-p pi0, omega -> pi+ pi- pi0
# Full final state: p anti-p pi+ pi- 3 pi0 -> p anti-p pi+ pi- 6 gamma
# ============================================================

decay_card_omega_mode = <<~DECAYCARD
  Decay psi(2S)
  1.000 Sigma+ anti-Sigma- omega PHSP;
  Enddecay
  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay
  Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
  Enddecay
  Decay omega
  1.000 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_SigmaSigma_omega"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_omega_mode
  config.cross_section   = :default
end

alg_omega = Algorithm.new("PsipSigmaSigmaOmega")
alg_omega.set_header(["PsipSigmaSigmaOmegaAlg/PsipSigmaSigmaOmega.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_omega = Selection.new
  # Charged track selection: p, anti-p (from Sigma decays), pi+, pi- (from omega)
  # Sigma daughter tracks use Vr < 2 cm due to hyperon lifetime;
  # omega daughter tracks use standard Vr < 1 cm.
  # Use Vr 2.0 to accommodate both — tighter cut on non-Sigma tracks
  # will be applied in for_each below.
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        2.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  end
  # Photon selection: 3 pi0 -> 6 photons minimum
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=6"
  end
  # PID: identify protons and anti-protons first
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  end
  # Remove identified protons/anti-protons from charged track list
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  # Identify remaining tracks as pions
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  end
  # 7C kinematic fit: 4C (four-momentum) + 3 x 1C (pi0 mass constraints on 3 photon pairs)
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) do
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 45
  end
  # Competing hypothesis: 5 gamma (wrong photon count) — stores chi2 for ROOT-level veto
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma]) do
    constrain_four_momentum
  end
  # Competing hypothesis: 7 gamma (wrong photon count)
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) do
    constrain_four_momentum
  end

alg_omega
  .note(:sigma_vr_cut, "Sigma+/Sigma- daughter tracks require Vr < 2 cm (hypperon lifetime); omega daughter tracks require Vr < 1 cm. select_track uses Vr 2.0 for all; tighter cut deferred to ROOT")
  .note(:sigma_selection, "Sigma+Sigma- candidates selected by minimizing delta_omega = sqrt((M_p_pi0 - m_Sigma)^2 + (M_anti_p_pi0 - m_anti_Sigma)^2 + (M_pi+pi-pi0 - m_omega)^2); applied in ROOT using fitted four-momenta")
  .note(:jpsi_veto, "J/psi-related background vetoed: RM(pi+pi-) outside J/psi mass window; RM(pi0pi0) outside J/psi mass window")
  .note(:eta_veto, "eta-related background vetoed: M(pi0 pi0 pi0) outside eta mass window")
  .note(:photon_count_veto, "chi2_4c(6gamma) < chi2_4c(5gamma) AND chi2_4c(6gamma) < chi2_4c(7gamma) required to suppress wrong-photon-count backgrounds")
  .note(:sigma_mass_window, "Sigma+Sigma- 2D signal region: both M(p pi0) in [1176, 1197] MeV/c^2; sideband regions used for background subtraction")
  .note(:bkg_subtraction, "non-Sigma+Sigma- peaking background from 2D sidebands subtracted with normalization factor from sideband fits")
  .with_decay_card(decay_card_omega_mode)
  .apply(sel_omega)

alg_omega.execute_on([psip_data, psip_incMC, exMC_omega])

# ============================================================
# Mode II: psi(3686) -> Sigma+ anti-Sigma- phi
# Sigma+ -> p pi0, anti-Sigma- -> anti-p pi0, phi -> K+ K-
# Partial reconstruction: one pi0 reconstructed, the other treated as missing
# ============================================================

decay_card_phi_mode = <<~DECAYCARD
  Decay psi(2S)
  1.000 Sigma+ anti-Sigma- phi PHSP;
  Enddecay
  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay
  Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_SigmaSigma_phi"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_phi_mode
  config.cross_section   = :default
end

alg_phi = Algorithm.new("PsipSigmaSigmaPhi")
alg_phi.set_header(["PsipSigmaSigmaPhiAlg/PsipSigmaSigmaPhi.h"])
        .set_constant({"ECMS" => [:double, 3.686]})

sel_phi = Selection.new
  # Charged track selection: p, anti-p (from Sigma decays), K+, K- (from phi)
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        2.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  end
  # Photon selection: partial reconstruction needs >= 2 photons for one pi0
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  # PID: identify protons and kaons
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==1"
  end
  # Remove identified protons/anti-protons and kaons from charged track list
  .remove([:prp <= :chrgp, :prm <= :chrgn, :kp <= :chrgp, :km <= :chrgn])
  # Reconstruct one pi0 -> gamma gamma (1-C Kalman fit)
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  end
  # 2C kinematic fit with one pi0 missing:
  # 1C from reconstructed pi0 mass constraint (applied in kalman fit above)
  # 1C from missing pi0 mass constraint (miss_track_of below)
  .kinematic_fit([:prp, :prm, :kp, :km, :pi0]) do
    nominal
    miss_track_of(:pi0)
    chi2_cut 20
  end

alg_phi
  .note(:sigma_vr_cut, "Sigma+/Sigma- daughter tracks require Vr < 2 cm (hyperon lifetime); phi daughter tracks prefer Vr < 1 cm. select_track uses Vr 2.0 for all")
  .note(:sigma_selection, "Sigma+Sigma- candidates selected by minimizing delta_phi = sqrt((M_p_pi0 - m_Sigma)^2 + (M_anti_p_pi0 - m_anti_Sigma)^2); applied in ROOT using fitted four-momenta")
  .note(:jpsi_veto, "J/psi-related background vetoed: RM(pi0pi0) outside J/psi mass window")
  .note(:sigma_mass_window, "Sigma+Sigma- 2D signal region: both M(p pi0) in [1176, 1200] MeV/c^2; sideband regions used for background subtraction")
  .note(:bkg_subtraction, "non-Sigma+Sigma- peaking background from 2D sidebands subtracted; Delta-related background negligible; chi_c1,2 -> Lambda Lambda phi peaking background found negligible")
  .note(:partial_reco, "partial reconstruction mode: only one pi0 reconstructed; the second pi0 treated as missing particle in the 2C kinematic fit")
  .with_decay_card(decay_card_phi_mode)
  .apply(sel_phi)

alg_phi.execute_on([psip_data, psip_incMC, exMC_phi])