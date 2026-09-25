# ============================================================
# Paper: arXiv:2410.13478v1
# First observation of chi_c0 → Sigma+ Sigma- eta using psi(3686) radiative decays
# Also evidence for chi_c1 → Sigma+ Sigma- eta and chi_c2 → Sigma+ Sigma- eta
# Sigma+ → p pi0, Sigma- → pbar pi0, eta → gamma gamma
# Data: (27.12±0.14)×10^8 psi(3686) events (BOSS 709)
# ============================================================

### ============================================================
### Datasets — psi(3686)
### ============================================================

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

### ============================================================
### Decay Card — chi_c0,1,2 share the same final state
### psi(3686) → gamma chi_cJ, chi_cJ → Sigma+ Sigma- eta
### Sigma+ → p pi0, Sigma- → anti-p- pi0, eta → gamma gamma
### ============================================================

decay_chi_cJ_SSeta = <<~DECAYCARD
  Decay psi(3686)
  1.000 gamma chi_c0 PHSP;
  Enddecay
  Decay chi_c0
  1.000 Sigma+ anti-Sigma- eta PHSP;
  Enddecay
  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay
  Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
  Enddecay
  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

### ============================================================
### Exclusive MC
### ============================================================

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_chi_cJ_SSeta"
  config.related_dataset = psip_data
  config.events          = 3000000
  config.decay_card      = decay_chi_cJ_SSeta
  config.cross_section   = :default
end

### ============================================================
### Algorithm: psi(3686) → gamma chi_cJ → gamma Sigma+ Sigma- eta
### Final state: radiative photon + p + pbar + pi0 + pi0 + eta
### = 1 radiative gamma + 2 charged tracks + 6 photons from pi0,eta
### Minimum: 7 photons, >=1 positive track, >=1 negative track
### 6C kinematic fit: 4C + 2 pi0 mass constraints
### Sigma+ decay length: Vz < 15 cm, Vr < 2 cm for p/pbar
### ============================================================

alg_chi_cJ_SSeta = Algorithm.new("ChiCJSSEta")
alg_chi_cJ_SSeta.set_header(["ChiCJSSEtaAlg/ChiCJSSEta.h"])
  .set_constant({"ECMS" => [:double, 3.686]})
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 6C kinematic fit")
  .note(:radiative_photon_selection, "among all remaining photons after pi0 and eta reconstruction, the one giving the largest DeltaM value (see Eq.1) is assigned as the radiative photon from psi(3686)")
  .note(:sigma_mass_window, "Sigma+/- signal mass window: M(ppi0) in [1.174, 1.204] GeV/c^2; sidebands [1.139,1.169] and [1.209,1.239] GeV/c^2")
  .note(:eta_mass_window, "eta signal mass window: M(gamma gamma) in [0.517, 0.577] GeV/c^2; sidebands [0.448,0.508] and [0.588,0.648] GeV/c^2")
  .note(:jpsi_veto, "veto chi_cJ → gamma J/psi: |M(Sigma+Sigma-gamma_E) - 3.091| > 0.056 GeV/c^2")
  .note(:pi0_veto, "veto chi_cJ → Sigma+Sigma-pi0: |M(gamma_E gamma1) - 0.132| > 0.018 and |M(gamma_E gamma2) - 0.137| > 0.044 GeV/c^2")
  .note(:background_veto, "chi2_signal < chi2_bkg for background channels formed by adding or subtracting one photon; competing 4C hypothesis fits for background rejection")
  .note(:signal_yield_fit, "simultaneous fit to M(Sigma+Sigma-gamma gamma) in eta signal and sideband regions; chi_c0 significance 7.0 sigma, chi_c1 4.3 sigma, chi_c2 4.6 sigma")

sel_chi_cJ_SSeta = Selection.new
sel_chi_cJ_SSeta.select_track {
  cos_theta 0.93
  Vz 15.0
  Vr 2.0
  nChrp ">=1"
  nChrn ">=1"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=7"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp ">=1"
  nprm ">=1"
}
.remove([:prp <= :chrgp, :prm <= :chrgn])
.assign({chrgp: :pip, chrgn: :pim})
# Reconstruct pi0 candidates from photon pairs — at least 2 needed (one for Sigma+ p pi0, one for Sigma- pbar pi0)
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0 ">=2"
}
# Main 6C kinematic fit (4C + 2 pi0 mass constraints): gamma_rad + p + pbar + 2 pi0 + eta
# Note: the radiative photon and eta photon pair selection (DeltaM method) are post-fit procedures — captured as notes
.kinematic_fit([:gamma, :prp, :prm, :pi0, :pi0, :gamma, :gamma]) {
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 45
}
# Competing-hypothesis background veto: 4C fit with extra/missing photon
# Signal hypothesis chi2_signal must be < chi2_bkg for all competing bg hypotheses
# This is captured as a note since the detailed competing hypothesis enumeration is post-fit ROOT logic
.kinematic_fit([:gamma, :prp, :prm, :pi0, :pi0, :gamma, :gamma, :gamma]) {
  constrain_four_momentum
}

alg_chi_cJ_SSeta.with_decay_card(decay_chi_cJ_SSeta).apply(sel_chi_cJ_SSeta)
alg_chi_cJ_SSeta.execute_on([psip_data, psip_incMC, exMC_signal])