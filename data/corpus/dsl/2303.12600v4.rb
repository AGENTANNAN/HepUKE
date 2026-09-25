# Paper: 2303.12600v4
# Title: Search for D_s+ -> tau+ nu_tau via tau+ -> pi+ nu_tau_bar with BDT method
# Energy: 8 energy points: 4.128-4.226 GeV
# e+e- -> D_s*+ D_s- pair production
# Double-tag technique: 13 hadronic ST modes for D_s-, signal D_s+ -> tau+ nu, tau+ -> pi+ nu_bar

### Dataset preparation ###
# 8 energy points at D_s* threshold: 4.128, 4.157, 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV
data_703_4128 = DatasetManager.load_real_data.find("703_4128")
incMC_703_4128 = DatasetManager.load_inclusive_mc.find("703_4128")

all_data = [data_703_4128]
all_incMC = [incMC_703_4128]

# Decay card: D_s*+ D_s- pair production via psi(4260) as KKMC top mother
# D_s*+ -> gamma D_s+; D_s+ -> tau+ nu_tau; tau+ -> pi+ nu_tau_bar
# D_s- decays inclusively (tag side, not specified in signal MC)
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000  D_s*+  D_s-                          PHSP;
    Enddecay

    Decay D_s*+
    1.000  gamma  D_s+                          VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.000  tau+  nu_tau                         PHSP;
    Enddecay

    Decay tau+
    1.000  pi+  anti-nu_tau                     PHSP;
    Enddecay

    Decay D_s-
    1.000  K+  K-  pi-                          PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Ds_tau_nu_pi_signal"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
# TagAnalysis: double-tag technique with D_s- hadronic ST and D_s+ -> tau+ nu signal
alg = TagAnalysis.new("DsTauNuBDT")
alg.set_header(["DsTauNuAlg/DsTauNuBDT.h"])
   .set_constant({"ECMS" => [:double, 4.178]})

# D_s- hadronic tag side (13 modes from paper Table 2)
# Known modes from DTagAlg; unavailable modes noted below
ds_tag_modes = [
  :DstoKKPi,        # K+ K- pi-
  :DstoKKPiPi0,     # K+ K- pi- pi0
  :DstoPiPiPi,      # pi+ pi- pi-
  :DstoKsK,         # K_S0 K-
  :DstoKsKPiPi,     # K_S0 K+ pi- pi-
  :DstoPiEta,       # pi- eta (eta -> gamma gamma)
  :DstoPiPi0Eta,    # rho- eta (rho- -> pi- pi0)
  :DstoEtaPPiPiEta, # pi- eta' (eta' -> pi+ pi- eta)
  :DstoEtaPRhoGam   # pi- eta' (eta' -> gamma rho0)
]

alg.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.charm -1
end

# Signal side: transition photon + pi+ from tau+ -> pi+ nu_bar
# Missing: combined neutrino system (nu_tau + nu_tau_bar)
alg.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 1   # opposite to D_s- tag (charm -1)
  s.missing :nu        # massless neutrino system
end

# 4C kinematic fit: tag D_s- + transition photon + pi+ + missing = ecms_lab
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.with_decay_card(decay_card).apply

# Tag modes unavailable in DTagAlg v1 — see tag-mode-vocabulary-gap memory
alg.note(:tag_mode_unavailable,
  "Four tag modes unavailable in DTagAlg v1: K_S0 K- pi0, K- pi+ pi-, K- K+ pi+ pi-, pi- eta(eta->3pi). Dropped from tag_side; these modes contribute ~20% of total ST yield. Applied in ROOT.")

# Transition photon from D_s*+ -> D_s+ gamma:
# - DeltaE = E_cm - E_tag - E_miss - E_gamma(pi0) minimized
# - Transition photon energy in D_s* rest frame: 0.114 < E_gamma < 0.149 GeV
# - Direct/indirect tag hypothesis: closest to nominal D_s* mass
alg.note(:transition_photon,
  "Transition photon from D_s*+ -> D_s+ gamma: |DeltaE| minimized; E_gamma in D_s* rest frame in [0.114,0.149] GeV; direct/indirect tag hypothesis selected by closest D_s* mass. Photon efficiency ~85%. Applied in ROOT.")

# Recoil mass M_rec against tagged D_s-: energy-dependent requirements (Table 1)
# Retain D_s*+ D_s- events; choose candidate with M_rec closest to D_s* mass
alg.note(:recoil_mass,
  "M_rec against tagged D_s- with E_cm-dependent requirements (Table 1). Candidate with M_rec closest to nominal D_s* mass retained. Applied in ROOT.")

# DT signal selection:
# - N_extra_char = 0 (exactly one additional pi+ track)
# - N_extra_pi0 = 0 (no additional pi0)
# - Pion PID for signal track
# - E/p < 0.9 (electron suppression)
# - Max extra photon energy E_neu_max < 0.3 GeV
# - |cos(theta_miss)| < 0.9 (point to fiducial volume)
# - M_miss^2 in [-0.2, 0.6] GeV^2/c^4
alg.note(:dt_signal_selection,
  "DT selection: N_extra_char=0; N_extra_pi0=0; pion PID; E/p<0.9; E_neu_max<0.3 GeV; |cos(theta_miss)|<0.9; M_miss^2 in [-0.2,0.6] GeV^2/c^4. Applied in ROOT.")

# BDT multivariate analysis with 9 input variables
alg.note(:bdt_analysis,
  "BDT trained with TMVA: input variables M_miss^2, M_tag, m_BC_tag, cos(theta_miss_gamma), cos(theta_miss), E_gamma(pi0), E_gamma_sum, cos(theta_pi+), p_pi+. Hyperparameters optimized. Applied in ROOT.")

# Background composition:
# - D_s+ -> mu+ nu_mu (38.8%), other tau decays (15.3%)
# - qqbar continuum (9.3%), tau+tau- (4.0%)
# - D_s+ -> eta pi+ (2.3%), D_s+ -> K_L0 K+ (2.7%), D_s+ -> K0 pi+ (4.2%)
# - Mixed open-charm (~23%)
# Four control regions for background validation
alg.note(:background_modeling,
  "Background: D_s+->mu+nu (38.8%), other tau (15.3%), qqbar (9.3%), tau+tau- (4.0%), eta pi+ (2.3%), K_L0 K+ (2.7%), K0 pi+ (4.2%), mixed open-charm (~23%). Four control regions: mu_nu, tau_other, qq_tau_tau, eta_pi. Applied in ROOT.")

# 8 energy points
alg.note(:energy_points,
  "8 c.m. energies: 4.128, 4.157, 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV. 7.33 fb^-1 total. ST yields and efficiencies in paper Tables 2,3.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)