# BESIII DSL: e+e- → π+π-π0 ηc and Zc(3900)± → ρ±ηc
# ArXiv: 1906.00831v3
# ηc reconstructed in 9 hadronic decay modes
# 5 energy points: 4.226, 4.258, 4.358, 4.416, 4.600 GeV

### Dataset preparation ###

data_4226 = DatasetManager.real_data.find("703_4230")
data_4360 = DatasetManager.real_data.find("703_4360")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4600 = DatasetManager.real_data.find("703_4600")

data_points = [data_4226, data_4360, data_4420, data_4600]

# Decay card for e+e- → π+π-π0 ηc, ηc → p pbar
# Using psi(4260) as top mother (single energy: the main √s=4.226 GeV point)
decay_card_pipipi0_etac_ppbar = <<~DECAYCARD
  Decay psi(4260)
  1.000 pi+ pi- pi0 eta_c PHSP;
  Enddecay

  Decay eta_c
  1.000 p+ anti-p- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for signal (ppbar channel, for the 4.226 GeV main energy)
sig_pipipi0_etac_ppbar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_pipipi0_etac_ppbar"
  config.related_dataset = data_4226
  config.events          = 500_000
  config.decay_card      = decay_card_pipipi0_etac_ppbar
  config.cross_section   = :default
end

### Event selection (BOSS) — representative ppbar channel ###

alg_pipipi0_etac = Algorithm.new("PiPiPi0Etac")
alg_pipipi0_etac.set_header(["PiPiPi0EtacAlg/PiPiPi0Etac.h"])
                 .set_constant({ "ECMS" => [:double, 4.226] })

sel_pipipi0_etac = Selection.new
  # Select charged tracks: 2π+ + 2π- (or 1p + 1pbar + 1π+ + 1π-)
  .select_track {
    nChrp ">=1"
    nChrn ">=1"
    nTot ">=4"
    nNet "==0"
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
  }
  # Select photons: 2 from π0, plus any from ηc decay products
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  }
  # PID: identify pions and protons
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion, against: [:kaon, :proton]
    nprp ">=1"
    nprm ">=1"
    npip ">=1"
    npim ">=1"
  }
  # Remove identified from generic lists
  .remove([:prp <= :chrgp, :prm <= :chrgn, :pip <= :chrgp, :pim <= :chrgn])
  # 1C Kalman fit to reconstruct π0 → γγ
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # 4C kinematic fit: e+e- → π+π-π0 ηc
  .kinematic_fit([:pip, :pim, :prp, :prm, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_pipipi0_etac
  .note(:etac_nine_decay_modes, "η_c reconstructed in 9 hadronic modes: p pbar, 2(K+K-), K+K-π+π-, K+K-π0, p pbar π0, K_S^0 K±π∓, π+π-η, K+K-η, π+π-π0π0; nine modes fitted simultaneously with unbinned max likelihood")
  .note(:etac_signal_region, "η_c signal region M(hadrons) in [2.95, 3.02] GeV/c^2; sidebands [2.78, 2.92] and [3.05, 3.19] GeV/c^2")
  .note(:etac_breit_wigner_fit, "η_c signal: constant-width Breit-Wigner convolved with Crystal Ball function; relative yields among channels constrained by BFs; background: 2nd order Chebyshev polynomial")
  .note(:d_meson_veto, "D meson veto: D^0→K±π∓ (|M-M_D0|>24 MeV), D^0→K±π∓π0, D±→K±π∓π± (|M-M_D±|>10 MeV), D±→K_S^0π±, D±→K_S^0π±π0")
  .note(:continuum_veto, "K*(892) veto: |M(Kπ) - m(K*)| > 32 MeV; ω veto: |M(π+π-π0) - m(ω)| > 26 MeV; η veto: |M(π+π-π0) - m(η)| > 10 MeV")
  .note(:peaking_background, "peaking background from e+e- → π+π-hc estimated with 600k MC events; contribution subtracted from signal yield")
  .note(:zc_search, "Zc(3900)± → ρ± ηc: ρ signal region M(π±π0) in [0.675, 0.875] GeV/c^2; pion momentum < 0.8 GeV/c; Zc mass/width fixed to latest BESIII measurement")
  .note(:pion_mass_assignment, "track masses assigned as kaon/pion/proton depending on decay mode; best combination chosen by minimizing χ2_4C + χ2_1C + χ2_PID + χ2_vertex")
  .note(:isr_correction, "ISR correction (1+δ) assuming signal from Y(4260); vacuum polarization correction applied; cross-check with e+e- → π+π-J/ψ energy-dependent cross section")
  .note(:multi_energy_points, "Analysis also performed at √s = 4.258, 4.358, 4.416, 4.600 GeV; upper limits determined at 90% C.L. using Bayesian method")
  .with_decay_card(decay_card_pipipi0_etac_ppbar)
  .apply(sel_pipipi0_etac)

alg_pipipi0_etac.execute_on(data_points + [sig_pipipi0_etac])