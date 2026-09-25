# Paper: 2304.12159v1
# Title: First search for D_s*+ -> e+ nu_e
# Energy: 8 energy points: 4.128-4.226 GeV
# e+e- -> D_s- D_s*+ pair production
# Double-tag technique: 16 hadronic ST modes for D_s-, D_s*+ -> e+ nu_e semileptonic decay

### Dataset preparation ###
data_703_4128 = DatasetManager.load_real_data.find("703_4128")
incMC_703_4128 = DatasetManager.load_inclusive_mc.find("703_4128")

all_data = [data_703_4128]
all_incMC = [incMC_703_4128]

# Decay card: D_s*+ D_s- pair production via psi(4260) as KKMC top mother
# D_s*+ -> e+ nu_e; D_s- decays inclusively
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000  D_s*+  D_s-                          PHSP;
    Enddecay

    Decay D_s*+
    1.000  e+  nu_e                             PHSP;
    Enddecay

    Decay D_s-
    1.000  K+  K-  pi-                          PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "DsStar_enu_signal"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
# TagAnalysis: double-tag with D_s- hadronic ST and D_s*+ -> e+ nu_e signal
alg = TagAnalysis.new("DsStarENuDTag")
alg.set_header(["DsStarENuDTagAlg/DsStarENuDTag.h"])
   .set_constant({"ECMS" => [:double, 4.178]})

# D_s- hadronic tag side: 16 modes (paper lists)
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

# Signal side: e+ + missing nu_e (semileptonic)
# D_s*+ -> e+ nu_e; no transition photon needed (direct D_s*+ decay)
alg.signal_side do |s|
  s.charged(ep: 1)
  s.require_charge 1    # opposite to D_s- tag (charm -1)
  s.missing :nu_e        # massless neutrino
end

# 4C kinematic fit
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.with_decay_card(decay_card).apply

# Tag modes unavailable in DTagAlg v1
alg.note(:tag_mode_unavailable,
  "Seven tag modes unavailable in DTagAlg v1: K_S0 K- pi0, K- pi+ pi-, K_S0 K_S0 pi-, K_S0 K- pi+ pi-, eta(pi+pi-pi0) pi-, eta'(eta_gammagamma) pi-, eta'(gamma rho0) pi-, eta(gammagamma) rho-, eta(pi+pi-pi0) rho-, eta(gammagamma) pi+ pi- pi-. Dropped from tag_side. Applied in ROOT.")

# ST selection:
# - Charged tracks: |cos(theta)|<0.93, |Vxy|<1cm, |Vz|<10cm (20cm for K_S0 daughters)
# - Kaon/pion PID: L(K/pi) > L(pi/K) and L(K/pi) > L(e)
# - K_S0: pi+pi- with |M(pi+pi-)-M(K_S0)| < 12 MeV/c^2, decay length > 2*sigma
# - pi0/eta: gamma gamma with 115-150 (pi0) or 500-570 (eta) MeV/c^2; mass-constrained kinematic fit
# - rho: M(pi pi) in (570,970) MeV/c^2; eta': M(eta pi+pi-) in (946,970), M(gamma rho0) in (940,976)
# - Momentum of direct pions from D_s- > 100 MeV/c (D* suppression)
# - K_S0 mass veto for pi+pi-pi- and K-pi+pi- tag modes
alg.note(:st_selection,
  "ST D_s- from 16 hadronic modes. Keep candidate with M(D_s-) closest to PDG value. M(D_s-) within 3sigma of nominal. M_rec from ST D_s-: signal region optimized by S/sqrt(S+B). ST yields from fit to M_rec: MC shape conv. Gaussian + 2nd/3rd-order Chebyshev. Applied in ROOT.")

# DT signal selection:
# - Exactly one extra charged track identified as e+
# - e+ PID: L_e > 0.8*(L_e+L_pi+L_K) and L_e > 0.001
# - Bremsstrahlung recovery: add EMC showers within 10 deg of e+ direction
# - Signal extracted from M_miss^2 distribution
alg.note(:dt_signal_selection,
  "DT signal: exactly one extra e+ track. e+ PID: L_e>0.8*(L_e+L_pi+L_K), L_e>0.001. Bremsstrahlung recovery. M_miss^2 signal extraction: 2.9sigma significance. B(D_s*+->e+nu_e) = (2.1+1.2-0.9_stat +/- 0.2_syst)x10^-5. Applied in ROOT.")

# D_s*+ decay constant from measured BF and total width
alg.note(:decay_constant,
  "Decay constant f_D_s*+ from measured BF: (213.6+61.0-45.8_stat +/- 43.9_syst) MeV. Applied in ROOT.")

# 8 energy points
alg.note(:energy_points,
  "8 c.m. energies: 4.128, 4.157, 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV. 7.33 fb^-1 total.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)