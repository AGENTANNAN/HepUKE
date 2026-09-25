# DSL for arXiv:2411.19642v2
# Measurement of inclusive cross sections of prompt J/psi and psi(3686) production
# in e+e- annihilation from sqrt(s)=3.808 to 4.951 GeV (49 energy points, 22 fb^-1)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# =============================================================================
# Datasets — 49 energy points across BOSS 703, 705, 706, 707
# =============================================================================

# BOSS 703 scan points (3.808–4.600 GeV)
pts_703 = [
  DatasetManager.real_data.find("703_3810"),  # 3.808 GeV, 50.54 pb^-1
  DatasetManager.real_data.find("703_3900"),  # 3.896 GeV, 52.61 pb^-1
  DatasetManager.real_data.find("703_4009"),  # 4.008 GeV, 482.0 pb^-1
  DatasetManager.real_data.find("703_4090"),  # 4.085 GeV, 52.86 pb^-1
  DatasetManager.real_data.find("703_4180"),  # 4.178 GeV, 3189.0 pb^-1
  DatasetManager.real_data.find("703_4190"),  # 4.189 GeV, 526.7 + 43.33 pb^-1
  DatasetManager.real_data.find("703_4200"),  # 4.199 GeV, 526.0 pb^-1
  DatasetManager.real_data.find("703_4210"),  # 4.209 GeV, 517.1 + 54.95 pb^-1
  DatasetManager.real_data.find("703_4220"),  # 4.219 GeV, 514.6 + 54.60 pb^-1
  DatasetManager.real_data.find("703_4230"),  # 4.226 GeV, 1056.4 + 44.54 pb^-1
  DatasetManager.real_data.find("703_4237"),  # 4.236 GeV, 530.3 pb^-1
  DatasetManager.real_data.find("703_4245"),  # 4.242 GeV, 55.88 pb^-1
  DatasetManager.real_data.find("703_4246"),  # 4.244 GeV, 538.1 pb^-1
  DatasetManager.real_data.find("703_4260"),  # 4.258 GeV, 828.4 pb^-1
  DatasetManager.real_data.find("703_4270"),  # 4.267 GeV, 531.1 pb^-1
  DatasetManager.real_data.find("703_4280"),  # 4.278 GeV, 175.7 pb^-1
  DatasetManager.real_data.find("703_4310"),  # 4.308 GeV, 45.08 pb^-1
  DatasetManager.real_data.find("703_4360"),  # 4.358 GeV, 543.9 pb^-1
  DatasetManager.real_data.find("703_4390"),  # 4.387 GeV, 55.57 pb^-1
  DatasetManager.real_data.find("703_4420"),  # 4.416 GeV, 1043.9 + 46.80 pb^-1
  DatasetManager.real_data.find("703_4470"),  # 4.467 GeV, 111.09 pb^-1
  DatasetManager.real_data.find("703_4530"),  # 4.527 GeV, 112.12 pb^-1
  DatasetManager.real_data.find("703_4575"),  # 4.575 GeV, 48.93 pb^-1
  DatasetManager.real_data.find("703_4600"),  # 4.600 GeV, 586.9 pb^-1
]

# BOSS 705 scan points (4.129–4.397 GeV)
pts_705 = [
  DatasetManager.real_data.find("705_4130"),  # 4.129 GeV, 401.5 pb^-1
  DatasetManager.real_data.find("705_4160"),  # 4.158 GeV, 408.7 pb^-1
  DatasetManager.real_data.find("705_4290"),  # 4.288 GeV, 502.4 pb^-1
  DatasetManager.real_data.find("705_4315"),  # 4.313 GeV, 501.2 pb^-1
  DatasetManager.real_data.find("705_4340"),  # 4.338 GeV, 505.0 pb^-1
  DatasetManager.real_data.find("705_4380"),  # 4.378 GeV, 522.7 pb^-1
  DatasetManager.real_data.find("705_4400"),  # 4.397 GeV, 507.8 pb^-1
  DatasetManager.real_data.find("705_4440"),  # 4.437 GeV, 569.9 pb^-1
]

# BOSS 706 scan points (4.612–4.699 GeV)
pts_706 = [
  DatasetManager.real_data.find("706_4610"),  # 4.612 GeV, 103.65 pb^-1
  DatasetManager.real_data.find("706_4620"),  # 4.628 GeV, 521.53 pb^-1
  DatasetManager.real_data.find("706_4640"),  # 4.641 GeV, 551.65 pb^-1
  DatasetManager.real_data.find("706_4660"),  # 4.661 GeV, 529.43 pb^-1
  DatasetManager.real_data.find("706_4680"),  # 4.682 GeV, 1667.39 pb^-1
  DatasetManager.real_data.find("706_4700"),  # 4.699 GeV, 535.54 pb^-1
]

# BOSS 707 scan points (4.740–4.951 GeV)
pts_707 = [
  DatasetManager.real_data.find("707_4740"),  # 4.740 GeV, 163.87 pb^-1
  DatasetManager.real_data.find("707_4750"),  # 4.750 GeV, 366.55 pb^-1
  DatasetManager.real_data.find("707_4780"),  # 4.781 GeV, 511.47 pb^-1
  DatasetManager.real_data.find("707_4840"),  # 4.843 GeV, 525.16 pb^-1
  DatasetManager.real_data.find("707_4914"),  # 4.918 GeV, 207.82 pb^-1
  DatasetManager.real_data.find("707_4946"),  # 4.951 GeV, 159.28 pb^-1
]

all_scan_points = pts_703 + pts_705 + pts_706 + pts_707

# Corresponding inclusive MC datasets
all_incMC = all_scan_points.map { |pt|
  DatasetManager.inclusive_mc.find("#{pt.boss}_#{pt.sample_name}")
}

# =============================================================================
# Decay cards — KKMC + psi(4260) top mother for inclusive charmonium production
# =============================================================================

# Decay card for inclusive J/psi production: e+e- → J/psi X, J/psi → mu+ mu-
decay_jpsi_inclusive = <<~DECAYCARD
  Decay psi(4260)
  1.0000 J/psi X        PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu-        VLL;
  Enddecay
  Decay X
  1.0000 pi+ pi-        PHSP;
  Enddecay
  End
DECAYCARD

# Decay card for ISR return to J/psi: e+e- → gamma_ISR J/psi, J/psi → mu+ mu-
decay_isr_jpsi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma J/psi    PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu-        VLL;
  Enddecay
  End
DECAYCARD

# Decay card for inclusive psi(3686) production: e+e- → psi(3686) X,
# psi(3686) → pi+ pi- J/psi, J/psi → mu+ mu-
decay_psip_inclusive_mu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 psi(2S) X      PHSP;
  Enddecay
  Decay psi(2S)
  1.0000 pi+ pi- J/psi  JPIPI;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu-        VLL;
  Enddecay
  Decay X
  1.0000 pi0 pi0        PHSP;
  Enddecay
  End
DECAYCARD

# Decay card for inclusive psi(3686) production: J/psi → e+ e-
decay_psip_inclusive_e = <<~DECAYCARD
  Decay psi(4260)
  1.0000 psi(2S) X      PHSP;
  Enddecay
  Decay psi(2S)
  1.0000 pi+ pi- J/psi  JPIPI;
  Enddecay
  Decay J/psi
  1.0000 e+ e-          VLL;
  Enddecay
  Decay X
  1.0000 pi0 pi0        PHSP;
  Enddecay
  End
DECAYCARD

# Decay card for ISR return to psi(3686): e+e- → gamma_ISR psi(3686),
# psi(3686) → pi+ pi- J/psi, J/psi → mu+ mu-
decay_isr_psip_mu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma psi(2S)  PHSP;
  Enddecay
  Decay psi(2S)
  1.0000 pi+ pi- J/psi  JPIPI;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu-        VLL;
  Enddecay
  End
DECAYCARD

# Decay card for ISR return to psi(3686): J/psi → e+ e-
decay_isr_psip_e = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma psi(2S)  PHSP;
  Enddecay
  Decay psi(2S)
  1.0000 pi+ pi- J/psi  JPIPI;
  Enddecay
  Decay J/psi
  1.0000 e+ e-          VLL;
  Enddecay
  End
DECAYCARD

# Decay card for chi_cJ inclusive: e+e- → chi_cJ X, chi_cJ → gamma J/psi, J/psi → mu+ mu-
decay_chicj_inclusive = <<~DECAYCARD
  Decay psi(4260)
  1.0000 chi_c1 X       PHSP;
  Enddecay
  Decay chi_c1
  1.0000 gamma J/psi    PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu-        VLL;
  Enddecay
  Decay X
  1.0000 pi+ pi-        PHSP;
  Enddecay
  End
DECAYCARD

# =============================================================================
# Exclusive MC samples — signal MC for efficiency determination (KKMC)
# =============================================================================

# J/psi inclusive signal MC (for efficiency epsilon_J/psiX)
exMC_jpsi = DatasetManager.create_exclusive_mc_for(all_scan_points) do |config|
  config.sample_name   = "sig_jpsi_inclusive"
  config.events        = 100_000
  config.decay_card    = decay_jpsi_inclusive
  config.cross_section = :default
end

# ISR J/psi signal MC (for R_gammaISR_J/psi ratio and N_gammaISR_J/psi background)
exMC_isr_jpsi = DatasetManager.create_exclusive_mc_for(all_scan_points) do |config|
  config.sample_name   = "sig_isr_jpsi"
  config.events        = 100_000
  config.decay_card    = decay_isr_jpsi
  config.cross_section = :default
end

# psi(3686) inclusive signal MC — muon channel (for efficiency epsilon^mu_psiX)
exMC_psip_mu = DatasetManager.create_exclusive_mc_for(all_scan_points) do |config|
  config.sample_name   = "sig_psip_inclusive_mu"
  config.events        = 100_000
  config.decay_card    = decay_psip_inclusive_mu
  config.cross_section = :default
end

# psi(3686) inclusive signal MC — electron channel (for efficiency epsilon^e_psiX)
exMC_psip_e = DatasetManager.create_exclusive_mc_for(all_scan_points) do |config|
  config.sample_name   = "sig_psip_inclusive_e"
  config.events        = 100_000
  config.decay_card    = decay_psip_inclusive_e
  config.cross_section = :default
end

# ISR psi(3686) signal MC — muon channel (for R^mu_gammaISR_psi ratio)
exMC_isr_psip_mu = DatasetManager.create_exclusive_mc_for(all_scan_points) do |config|
  config.sample_name   = "sig_isr_psip_mu"
  config.events        = 100_000
  config.decay_card    = decay_isr_psip_mu
  config.cross_section = :default
end

# ISR psi(3686) signal MC — electron channel (for R^e_gammaISR_psi ratio)
exMC_isr_psip_e = DatasetManager.create_exclusive_mc_for(all_scan_points) do |config|
  config.sample_name   = "sig_isr_psip_e"
  config.events        = 100_000
  config.decay_card    = decay_isr_psip_e
  config.cross_section = :default
end

# chi_cJ inclusive signal MC (for efficiency epsilon_chi_cJ X)
exMC_chicj = DatasetManager.create_exclusive_mc_for(all_scan_points) do |config|
  config.sample_name   = "sig_chicj_inclusive"
  config.events        = 100_000
  config.decay_card    = decay_chicj_inclusive
  config.cross_section = :default
end

# =============================================================================
# Algorithm 1 — J/psi inclusive: e+e- → J/psi X, J/psi → mu+ mu-
# (Also used for e+e- → gamma_ISR J/psi, J/psi → mu+ mu- — same final state)
# =============================================================================

alg_jpsi = Algorithm.new("JpsiInclusive")
alg_jpsi.set_header(["JpsiInclusiveAlg/JpsiInclusive.h"])
alg_jpsi.note(:inclusive_yield_extraction, "signal yield N_J/psiX extracted by fitting di-muon
  invariant mass M(mu+mu-) distribution in range (2.8,3.4) GeV/c^2 with quadratic background
  function; signal region (3.0,3.2) GeV/c^2 excluded from background fit")
alg_jpsi.note(:isr_subtraction, "ISR return to J/psi background subtracted via observed
  e+e-→gamma_ISR J/psi events; ratio R_gammaISR_J/psi calculated from KKMC samples using Eq.(14)")
alg_jpsi.note(:feeddown_subtraction, "psi(3686)→J/psi X and chi_cJ→gamma J/psi feed-down
  contributions subtracted using Eqs.(3)-(6); other charmonium decays (psi(3770), chi_c0)
  neglected due to negligible contribution")
alg_jpsi.note(:born_cross_section, "Born cross section obtained via iterative ISR correction
  procedure Eq.(8)-(13) with vacuum polarization factor |1-Pi|^{-2}")

sel_jpsi = Selection.new
sel_jpsi.select_track {
  cos_theta 0.93     # |cos(theta)| < 0.93
  Vz 10.0            # |Vz| < 10 cm
  Vr 1.0             # |Vxy| < 1.0 cm
  nChrp ">=1"        # at least one positive track
  nChrn ">=1"        # at least one negative track
  nNet "==0"         # net charge zero
}
.select_photon {
  tdc_emc_start 0    # TDC start: 0 ns
  tdc_emc_end 14     # TDC end: 700 ns (14 x 50 ns)
  energyThreshold_b 0.025  # E > 25 MeV barrel
  energyThreshold_e 0.050  # E > 50 MeV endcap
  angle_to_track 20.0      # min angle to nearest charged track > 20 deg
  nGam ">=0"         # any number of photons (Table II)
}
.pid(method: :probability) {
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  nlp ">=1"          # at least one positive lepton
  nlm ">=1"          # at least one negative lepton
}
# Kalman fit: constrain di-muon invariant mass to nominal J/psi mass
# This improves invariant mass resolution for psi(3686) and chi_cJ reconstruction
.kalman_kinematic_fit([:mup, :mum]) {
  invariant_mass_of(:mup, :mum).constrain_to_nominal_mass_of(:jpsi)
  chi2_cut 50        # choose combination with minimum chi2
}
# Store fitted momenta and final 4C kinematic fit
.kinematic_fit([:mup, :mum]) {
  nominal
  constrain_four_momentum
  chi2_cut 200       # loose cut; tight cut in ROOT
}

alg_jpsi.with_decay_card(decay_jpsi_inclusive).apply(sel_jpsi)
alg_jpsi.execute_on(all_scan_points + all_incMC + exMC_jpsi + exMC_isr_jpsi)

# =============================================================================
# Algorithm 2 — psi(3686) inclusive: e+e- → psi(3686) X, psi(3686) → pi+ pi- J/psi,
#   J/psi → mu+ mu-
# =============================================================================

alg_psip_mu = Algorithm.new("PsipInclusiveMu")
alg_psip_mu.set_header(["PsipInclusiveMuAlg/PsipInclusiveMu.h"])
alg_psip_mu.note(:inclusive_yield_extraction, "signal yield N^mu_psiX extracted by fitting
  pi+pi-J/psi invariant mass M(pi+pi-mu+mu-) in range (3.50,3.90) GeV/c^2; signal region
  (3.65,3.72) GeV/c^2 excluded; double Gaussian + linear/quadratic background model
  depending on sqrt(s)")
alg_psip_mu.note(:isr_subtraction, "ISR return to psi(3686) background subtracted via observed
  e+e-→gamma_ISR psi(3686) events; ratio R^mu_gammaISR_psi from KKMC Eq.(15)")
alg_psip_mu.note(:exclusive_peak, "exclusive e+e-→pi+pi-psi(3686) peak at ~3.636 GeV/c^2
  characterized by double Gaussian with fixed width from phase space MC")

sel_psip_mu = Selection.new
sel_psip_mu.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"        # at least 2 positive tracks (pi+, mu+)
  nChrn ">=2"        # at least 2 negative tracks (pi-, mu-)
  nNet "==0"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 20.0
  nGam ">=0"         # any number of photons (Table II)
}
.pid(method: :probability) {
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  nlp ">=1"
  nlm ">=1"
}
# Remove identified leptons from charged lists; assign remaining as pions
.remove([:lp <= :chrgp])
.remove([:lm <= :chrgn])
.assign({ chrgp: :pip, chrgn: :pim })
# Kalman fit: constrain di-muon to J/psi mass
.kalman_kinematic_fit([:mup, :mum]) {
  invariant_mass_of(:mup, :mum).constrain_to_nominal_mass_of(:jpsi)
  chi2_cut 50
}
# 4C kinematic fit for psi(3686) → pi+ pi- J/psi, J/psi → mu+ mu-
.kinematic_fit([:pip, :pim, :mup, :mum]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_psip_mu.with_decay_card(decay_psip_inclusive_mu).apply(sel_psip_mu)
alg_psip_mu.execute_on(all_scan_points + all_incMC + exMC_psip_mu + exMC_isr_psip_mu)

# =============================================================================
# Algorithm 3 — psi(3686) inclusive: e+e- → psi(3686) X, psi(3686) → pi+ pi- J/psi,
#   J/psi → e+ e-
# =============================================================================

alg_psip_e = Algorithm.new("PsipInclusiveE")
alg_psip_e.set_header(["PsipInclusiveEAlg/PsipInclusiveE.h"])
alg_psip_e.note(:inclusive_yield_extraction, "signal yield N^e_psiX extracted by fitting
  pi+pi-J/psi invariant mass in di-electron channel; same fitting method as muon channel")
alg_psip_e.note(:lepton_averaging, "observed cross sections from di-muon (sigma^O_psi,mu) and
  di-electron (sigma^O_psi,e) modes averaged taking only independent uncertainties")

sel_psip_e = Selection.new
sel_psip_e.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
  nNet "==0"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 20.0
  nGam ">=0"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  nlp ">=1"
  nlm ">=1"
}
.remove([:lp <= :chrgp])
.remove([:lm <= :chrgn])
.assign({ chrgp: :pip, chrgn: :pim })
.kalman_kinematic_fit([:ep, :em]) {
  invariant_mass_of(:ep, :em).constrain_to_nominal_mass_of(:jpsi)
  chi2_cut 50
}
.kinematic_fit([:pip, :pim, :ep, :em]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_psip_e.with_decay_card(decay_psip_inclusive_e).apply(sel_psip_e)
alg_psip_e.execute_on(all_scan_points + all_incMC + exMC_psip_e + exMC_isr_psip_e)

# =============================================================================
# Algorithm 4 — chi_cJ inclusive: e+e- → chi_cJ X, chi_cJ → gamma J/psi, J/psi → mu+ mu-
# =============================================================================

alg_chicj = Algorithm.new("ChicJInclusive")
alg_chicj.set_header(["ChicJInclusiveAlg/ChicJInclusive.h"])
alg_chicj.note(:inclusive_yield_extraction, "signal yields N_chi_c1X and N_chi_c2X extracted
  by fitting gamma J/psi invariant mass in range (3.45,3.70) GeV/c^2; exponential + Gaussian
  background model; signal regions (3.47,3.53) and (3.53,3.60) GeV/c^2 for chi_c1 and chi_c2
  respectively excluded from fit")
alg_chicj.note(:isr_psip_peak, "for sqrt(s) in [4.085,4.397] GeV, Gaussian peak from
  e+e-→gamma_ISR psi(3686) with gamma J/psi including ISR photon accounted for in background;
  Gaussian mean and width from KKMC samples")
alg_chicj.note(:assumed_source, "chi_cJ production assumed dominated by ISR return to
  psi(3686) resonance; efficiency evaluated from KKMC e+e-→gamma_ISR psi(3686) with
  psi(3686)→gamma chi_cJ→gamma(gamma J/psi), J/psi→mu+mu-")

sel_chicj = Selection.new
sel_chicj.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
  nChrn ">=1"
  nNet "==0"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 20.0
  nGam ">=1"         # at least 1 photon for chi_cJ → gamma J/psi
}
.select_isolated_photon {
  angle_to_mum_track 20.0   # isolation from muon tracks
  angle_to_mup_track 20.0
  nGam ">=1"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  nlp ">=1"
  nlm ">=1"
}
.kalman_kinematic_fit([:mup, :mum]) {
  invariant_mass_of(:mup, :mum).constrain_to_nominal_mass_of(:jpsi)
  chi2_cut 50
}
# 4C kinematic fit: chi_cJ → gamma J/psi, J/psi → mu+ mu-
.kinematic_fit([:gamma, :mup, :mum]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}
# Competing hypothesis: ISR psi(3686) → gamma J/psi with ISR photon as the gamma
# (store chi2 for ROOT-level veto)
.kinematic_fit([:gamma, :mup, :mum]) {
  invariant_mass_of(:mup, :mum).constrain_to_nominal_mass_of(:jpsi)
  constrain_four_momentum
}

alg_chicj.with_decay_card(decay_chicj_inclusive).apply(sel_chicj)
alg_chicj.execute_on(all_scan_points + all_incMC + exMC_chicj)