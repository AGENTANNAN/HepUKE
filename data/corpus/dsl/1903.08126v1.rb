# BOSS DSL for: e+e- -> pi+pi-DDbar (D = D0 or D+) cross section measurement
# Paper: arXiv:1903.08126v1 (BESIII, multi-energy scan above 4.08 GeV)
# Studies: pi+pi-psi(3770), rho0 X2(4013), D1(2420)Dbar in the pi+pi-DDbar final state
# D mesons fully reconstructed in hadronic decay modes (NOT DTagAlg — ordinary analysis)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Energy points: 4.09, 4.19, 4.21, 4.22, 4.23, 4.245, 4.26, 4.31, 4.36, 4.39, 4.42, 4.47, 4.53, 4.575, 4.60 GeV
# Matching BOSS 703 samples: 4090, 4190, 4210, 4220, 4230, 4245/4246, 4260, 4310, 4360, 4390, 4420, 4470, 4530, 4575, 4600
# Use highest-luminosity entry where multiple rounds exist

scan_points = [
  DatasetManager.real_data.find("703_4090"),   # ~4.085 GeV
  DatasetManager.real_data.find("703_4190"),   # ~4.189 GeV (round10, highest lumi)
  DatasetManager.real_data.find("703_4210"),   # ~4.208/4.209 GeV (round10)
  DatasetManager.real_data.find("703_4220"),   # ~4.217/4.219 GeV (round10)
  DatasetManager.real_data.find("703_4230"),   # ~4.226 GeV (round06 main, 1056/pb)
  DatasetManager.real_data.find("703_4246"),   # ~4.244 GeV (round10, closest to 4.242)
  DatasetManager.real_data.find("703_4260"),   # ~4.258 GeV
  DatasetManager.real_data.find("703_4310"),   # ~4.308 GeV
  DatasetManager.real_data.find("703_4360"),   # ~4.358 GeV
  DatasetManager.real_data.find("703_4390"),   # ~4.387 GeV
  DatasetManager.real_data.find("703_4420"),   # ~4.416 GeV (round07, 1044/pb)
  DatasetManager.real_data.find("703_4470"),   # ~4.467 GeV
  DatasetManager.real_data.find("703_4530"),   # ~4.527 GeV
  DatasetManager.real_data.find("703_4575"),   # ~4.575 GeV
  DatasetManager.real_data.find("703_4600"),   # ~4.600 GeV
]

# ============================================================
# Signal channel 1: e+e- -> pi+ pi- D0 anti-D0
# D0 -> K- pi+, anti-D0 -> K+ pi- (simplest mode, used as primary)
# ============================================================

decay_card_d0_kpi = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC batch for all energy points
exMC_d0_kpi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name = "pipiD0D0bar_D0toKPi"
  config.events = 200_000
  config.decay_card = decay_card_d0_kpi
  config.cross_section = :default
end

# Algorithm 1: D0 -> K- pi+, anti-D0 -> K+ pi- (both simplest)
alg_d0_kpi = Algorithm.new("PipiD0D0barD0toKPi")
alg_d0_kpi.set_header(["PipiD0D0barD0toKPiAlg/PipiD0D0barD0toKPi.h"])
            .set_constant({ "ECMS" => [:double, 4.260] })

sel_d0_kpi = Selection.new

sel_d0_kpi.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"    # at least K+ and pi+ from D0bar + pi+ from ISR
  nChrn ">=2"    # at least K- and pi- from D0 + pi- from ISR
end

sel_d0_kpi.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 20.0    # paper: > 20 degrees from nearest charged track
  nGam ">=0"
end

sel_d0_kpi.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion]
  nkp ">=1"
  nkm ">=1"
end

sel_d0_kpi.remove([:kp <= :chrgp, :km <= :chrgn])

sel_d0_kpi.assign({ chrgp: :pip, chrgn: :pim })

# 4C kinematic fit: pi+ pi- K+ K- equal to total pi+pi-KK system
# In the full analysis, this includes all decay products; simplified here for the main mode
sel_d0_kpi.kinematic_fit([:pip, :pim, :kp, :km, :pip, :pim]) do
  nominal
  constrain_four_momentum
  chi2_cut 56     # paper: chi2 < 56 for D0D0bar mode (optimised by FOM)
end

alg_d0_kpi.note(:multi_energy, "Analysis performed at 15 energy points from 4.0854 to 4.5995 GeV. The ECMS constant is set to 4.260 GeV representative; actual per-run CMS energy read from MeasuredEcmsSvc at execute time.")
   .note(:d_decay_modes, "D0 reconstructed in 4 modes (K-pi+, K-pi+pi0, K-pi+pi+pi-, K-pi+pi+pi-pi0); D+ in 5 modes (K-pi+pi+, K-pi+pi+pi0, K_S0 pi+, K_S0 pi+pi0, K_S0 pi+pi-pi+). Only the simplest mode (D0->K-pi+) is expressed in this DSL. Other modes require separate Algorithm objects with different track/photon/PID requirements.")
   .note(:pi0_reconstruction, "pi0 -> gamma gamma with M(gammagamma) in [0.115, 0.150] GeV/c^2; 1C kinematic fit constraining to pi0 mass for modes with pi0. Not expressed in this simplified DSL for the K-pi+ only mode.")
   .note(:ks_reconstruction, "K_S0 -> pi+pi- with vertex fit chi2 < 100 and M(pipi) in [0.487, 0.511] GeV/c^2 for modes with K_S0. Not expressed in this simplified DSL.")
   .note(:d_meson_selection, "Best DDbar pair chosen by average mass closest to nominal D mass. Signal region: -6 < DeltaMhat < 10 MeV/c^2 and |DeltaM| < 35 MeV/c^2 for D0D0bar, -5 < DeltaMhat < 10 MeV/c^2 and |DeltaM| < 25 MeV/c^2 for D+D-. DeltaMhat = (M(D)+M(Dbar))/2 - m_D, DeltaM = M(D)-M(Dbar).")
   .note(:dstar_veto, "Veto D* signal: M(D0 pi+) > 2.017 GeV/c^2 and M(D0bar pi-) > 2.017 GeV/c^2 for pi+pi-psi(3770) and rho0 X2(4013) channels. For D1(2420)0 -> D*+pi-: M(D0 pi+) < 2.017 and M(D0bar pi-) > 2.017 (and cc).")
   .note(:signal_extraction, "Signal yields from unbinned maximum likelihood fit to M(DDbar) with non-parametric kernel-estimation signal shape (from MC), D1(2420)Dbar + 4-body pipiDDbar + non-DDbar sideband backgrounds. X2(4013) search uses 3rd-order polynomial background.")
   .note(:cross_section_formula, "Born cross section: sigmaB = Nobs / (L_int * f_r * f_v * (B_N * sum(e_ij*B_i*B_j) + B_C * sum(e_kl*B_k*B_l))) where f_r = (1+delta^r) ISR correction, f_v = 1/|1-Pi|^2 vacuum polarization factor, B_N/B_C are psi(3770)->D0D0bar/D+D- BFs. X2(4013) and D1(2420)Dbar cross sections use similar formulas.")
   .note(:radiative_correction, "ISR correction factor f_r computed iteratively from QED radiator function F(x,s). Y(4260) line shape assumed as initial input; iteration continued until difference between iterations comparable to statistical uncertainty.")

alg_d0_kpi.with_decay_card(decay_card_d0_kpi).apply(sel_d0_kpi)
alg_d0_kpi.execute_on(scan_points + exMC_d0_kpi)