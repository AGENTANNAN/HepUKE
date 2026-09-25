# Paper: 2303.12927v1
# Title: D_s+ -> f0(980) e+ nu_e with f0(980) -> pi+ pi-
# Energy: 8 energy points: 4.128-4.226 GeV
# Double-tag technique: 12 hadronic ST modes for D_s-, signal D_s+ -> f0(980) e+ nu_e
# Semileptonic decay with missing neutrino

### Dataset preparation ###
data_703_4128 = DatasetManager.load_real_data.find("703_4128")
incMC_703_4128 = DatasetManager.load_inclusive_mc.find("703_4128")

all_data = [data_703_4128]
all_incMC = [incMC_703_4128]

# Decay card: D_s*+ D_s- pair production via psi(4260) as KKMC top mother
# D_s*+ -> gamma D_s+; D_s+ -> f0(980) e+ nu_e; f0(980) -> pi+ pi-
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000  D_s*+  D_s-                          PHSP;
    Enddecay

    Decay D_s*+
    1.000  gamma  D_s+                          VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.000  f_0  e+  nu_e                        PHSP;
    Enddecay

    Decay f_0
    1.000  pi+  pi-                             PHSP;
    Enddecay

    Decay D_s-
    1.000  K+  K-  pi-                          PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Ds_f0_enu_signal"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
# TagAnalysis: double-tag technique with D_s- hadronic ST and D_s+ -> f0(980) e+ nu_e signal
alg = TagAnalysis.new("DsF0ENu")
alg.set_header(["DsF0ENuAlg/DsF0ENu.h"])
   .set_constant({"ECMS" => [:double, 4.178]})

# D_s- hadronic tag side: 12 modes (paper)
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

# Signal side: pi+ pi- e+ + transition photon + missing neutrino
# Lepton key :ep drives electron ID
alg.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)
  s.require_charge 0   # +1 -1 +1 = +1 for D_s+, balanced by missing nu
  s.missing :nu_e       # massless neutrino
end

# 4C kinematic fit + f0(980) mass constraint
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.with_decay_card(decay_card).apply

# Tag modes unavailable in DTagAlg v1
alg.note(:tag_mode_unavailable,
  "Three tag modes unavailable in DTagAlg v1: K_S0 K- pi0, K- pi+ pi-, K_S0 K- pi+ pi-. Dropped from tag_side. Applied in ROOT.")

# Tag selection details: same as Ref [37]
# M_rec against tag D_s- selects e+e- -> D_s*+ D_s- events
# Multiple candidates: keep one with M_rec closest to nominal D_s* mass
# M_tag fits: signal = MC shape conv. Gaussian; BG = Chebyshev polynomial
alg.note(:tag_selection,
  "Tag selection: same as Ref [37]. M_rec against tag D_s-; candidate with M_rec closest to D_s* mass kept. M_tag fit: MC shape conv. Gaussian + Chebyshev polynomial. Applied in ROOT.")

# Signal side: 3 tracks (pi+, pi-, e+) + transition photon reconstruction
# M_rec^2 against transition photon + tag D_s-: must be in (3.78, 4.05) GeV^2/c^4
# |M_miss^2| < 0.06 GeV^2/c^4
alg.note(:signal_selection,
  "Signal: pi+ pi- e+ with transition photon. M_rec^2 vs transition gamma + tag D_s- in [3.78,4.05] GeV^2/c^4. |M_miss^2| < 0.06 GeV^2/c^4. M(pi+pi-) in [0.6,1.6] GeV/c^2. Applied in ROOT.")

# f0(980) lineshape modeled by Flatté formula (BESII parameters)
# Weighted signal efficiency: (35.44 +/- 0.07)%
# Background from D_s+ -> eta'(gamma pi+pi-) e+ nu_e peaking around 0.75 GeV/c^2
alg.note(:f0_signal_extraction,
  "f0(980) signal from unbinned ML fit to M(pi+pi-). Signal: MC shape conv. Gaussian; BG: inclusive MC shape conv. Gaussian. 439+/-33 signal events. B(D_s+->f0(980)e+nu, f0->pi+pi-) = (1.72+/-0.13_stat+/-0.10_syst)x10^-3. Applied in ROOT.")

# q^2-dependent FF extraction using simple pole parameterization
# 4 q^2 intervals; efficiency matrix from signal MC
alg.note(:form_factor_extraction,
  "Form factor f_+^f0(q^2) from decay rate in 4 q^2 intervals. Simple pole parameterization (M_pole=2.46 GeV/c^2). f_+^f0(0)|V_cs| = 0.504+/-0.017_stat+/-0.035_syst. Applied in ROOT.")

# 8 energy points
alg.note(:energy_points,
  "8 c.m. energies: 4.128, 4.157, 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV. 7.33 fb^-1 total. Total tag yield: 771101+/-3445.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)