# BOSS DSL for: J/psi -> K+ K- pi0 Partial-Wave Analysis
# Paper: arXiv:1904.10630v2 (BESIII, 223.7e6 J/psi events from 2009)
# PWA reveals: K*(892)+-, K2*(1430)+-, K2*(1980)+-, K4*(2045)+- in Kpi channel
#              rho-like JPC=1-- structures at ~1.65 and ~2.05 GeV/c^2 in KK channel
# pi0 -> gamma gamma

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for J/psi -> K+ K- pi0
decay_card = <<~DECAYCARD
    Decay J/psi
    1.000 K+ K- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "JpsiToKKPi0_PWA"
  config.related_dataset = data_jpsi
  config.events = 1_000_000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("JpsiKKPi0_PWA")
alg.set_header(["JpsiKKPi0_PWAAlg/JpsiKKPi0_PWA.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

sel = Selection.new

# Charged track selection: 2 tracks, net charge zero
# Paper: |cos theta| < 0.93, |z| < 10 cm, R < 1 cm, p_T > 120 MeV/c
sel.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
  nChrn ">=1"
  nNet "==0"
end

# Photon selection: at least 2 photons for pi0 -> gamma gamma
# Paper: E > 25 MeV barrel, E > 50 MeV endcap, angle > 10 deg, 0 < t < 700 ns
sel.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=2"
end

# PID: both charged tracks identified as kaons
# Paper: combined dE/dx + TOF C.L. for pi, K, p; highest C.L. assigned
sel.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==1"
  nkm "==1"
end

sel.remove([:kp <= :chrgp, :km <= :chrgn])

# 5C kinematic fit: 4-momentum conservation + pi0 mass constraint
# This is the final fit used for PWA momenta
# Paper: chi2_4C(gamma gamma K+ K-) < 60; then 5C fit constrains M(gg) to m_pi0
sel.kinematic_fit([:kp, :km, :gamma, :gamma]) do
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 60
  nominal
end

# Inexpressible BOSS-side procedures captured as notes
alg.note(:transverse_momentum_cut, "Events rejected if any charged track has p_T < 120 MeV/c. The DSL select_track does not support per-track p_T cuts; applied at ROOT level.")
   .note(:pid_highest_cl, "PID uses combined dE/dx + TOF confidence levels for pi, K, p hypotheses. Track assigned particle type with highest C.L. Both tracks must be identified as kaons. The DSL probability PID with identify :kaon approximates this selection.")
   .note(:photon_timing, "EMC cluster timing 0 <= t <= 700 ns from event start. DSL tdc_emc_start 0 / tdc_emc_end 14 (14 x 50 ns = 700 ns) approximates this.")
   .note(:diphoton_preselection, "Before the 4C fit, only photon pairs with M(gamma gamma) < 300 MeV/c^2 are considered. The DSL kinematic_fit automatically iterates all gamma-gamma combinations; the M_gg < 300 MeV/c^2 preselection is applied at ROOT level.")
   .note(:background_hypothesis_rejection, "After the gamma gamma K+ K- 4C fit, event compared against background hypotheses: gamma gamma pi+ pi- (pion mis-ID), gamma K+ K- (missing photon), gamma gamma gamma K+ K- (extra photon). If any background hypothesis gives lower chi2, event rejected. This multi-hypothesis comparison is performed at ROOT level.")
   .note(:pi0_mass_window, "pi0 candidates required: 110 < M(gamma gamma) < 150 MeV/c^2 after 4C preselection. Applied at ROOT level.")
   .note(:fivec_fit, "The 5C kinematic fit (4C + pi0 mass constraint) provides final particle momenta used in PWA. The DSL expresses this via constrain_four_momentum + invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0) inside the kinematic_fit block.")
   .note(:pwa_method, "Partial-wave analysis: isobar model with covariant tensor formalism. Event-by-event unbinned maximum likelihood fit over 5D phase space. Blatt-Weisskopf barrier factors at each vertex. 26+ intermediate amplitudes across K+ pi0, K- pi0, and K+ K- channels. Detector resolution correction via approximation for K*(892). Two solutions reported (I: well-established resonances only; II: with non-resonant J^P=3^- Kpi contribution). PWA and related ROOT-level fits are beyond the DSL scope.")
   .note(:resonances_observed, "Stable contributions: K*(892)+- (M~894, Gamma~47, ~89% fraction), K2*(1430)+- (M~1431, Gamma~100, ~9%), K2*(1980)+- (M~1817-1868, Gamma~272-312, ~0.4%), K4*(2045)+- (M~2015-2090, Gamma~183-201, ~0.16%). K+K- channel: JPC=1-- at ~1650 MeV/c^2 (Gamma~167-194, possible rho(1700)/omega(1650)/interference) and at ~2050 MeV/c^2 (Gamma~149-193, possible rho(2150)). Marginal: K*(1410)+-, K*(1680)+-, K3*(1780)+-, rho(770), rho(1450), rho3(1690).")
   .note(:branching_fraction, "B(J/psi -> K+ K- pi0) = (2.88 +/- 0.01(stat) +/- 0.12(syst)) x 10^-3. Selection efficiency from PWA solution II weighted MC. Branching fractions for individual resonant contributions also reported. Systematic 4.3% from: N_bg, N_continuum, tracking 2%, PID 2%, photon reco 2%, kinematic fit cut 2.4%, N_J/psi 0.6%.")
   .note(:background_estimation, "Total background 0.3% (565 +/- 24 from inclusive MC). Peaking background from J/psi -> gamma eta_c, eta_c -> K+ K- pi0. Non-peaking from pi0 mass sideband 190 < M_gg < 230 MeV/c^2. Continuum from 3.08 GeV data (~280 pb^-1): 855 +/- 499 events. Background subtracted from NLL in PWA.")
   .note(:dataset_info, "Uses (223.7 +/- 1.4) x 10^6 J/psi events from 2009 data-taking. DSL uses standard 708_3097 dataset. Also continuum data at 3.08 GeV (~280 pb^-1) for background estimation (separate dataset, not included in this DSL).")

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on([data_jpsi, incMC_jpsi, exMC_signal])