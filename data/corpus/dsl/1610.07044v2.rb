# =============================================================================
# 1610.07044v2
# Measurement of the cross section of e+e- -> pi+ pi- h_c at c.m. energies
# from 3.896 to 4.600 GeV (BESIII)
#
# The h_c is reconstructed through its electric-dipole transition
# h_c -> gamma eta_c, with eta_c -> X_i, where X_i is one of 16 exclusive
# hadronic final states: p pbar, 2(pi+ pi-), 2(K+ K-), pi+ pi- K+ K-,
# pi+ pi- p pbar, 3(pi+ pi-), 2(pi+ pi-) K+ K-, K_S0 K+- pi-+,
# K_S0 K+- pi-+ pi+ pi-, K+ K- pi0, p pbar pi0, K+ K- eta, pi+ pi- eta,
# 2(pi+ pi-) eta, pi+ pi- pi0 pi0 and 2(pi+ pi- pi0).
#
# BOSS part only : dataset preparation + event selection up to the final 4C
# kinematic fit.  The extraction of the number of signal events from the
# M(gamma eta_c) and M(pi+ pi- h_c) spectra, the Born cross sections and the
# coherent sum of two Breit-Wigner functions fitted to the dressed cross
# section (Y(4220), Y(4390)) are ROOT-level steps.
#
# Generator : the paper simulates ISR with KKMC (maximum ISR-photon energy
# set by the pi+ pi- h_c mass threshold) and generates e+e- -> pi+ pi- h_c
# according to phase space, so the standard KKMC + psi(4260) top-mother
# convention is used (not ConExc).
#
# The 16 eta_c decay modes have different final states and different selection
# chains, so (Rule T1) each mode gets its own Algorithm object and Selection.
# =============================================================================

# =============================================================================
# Datasets : the 17 c.m. energy points with an integrated luminosity larger
# than 40 pb^-1 ("XYZ data sample", Table II).  The corresponding sample names
# are taken from the BESIII dataset table (BOSS 703).
# The 62 "R-scan" points (10 MeV steps between 4.097 and 4.587 GeV,
# luminosities below 20 pb^-1, Table III) are stored in the same 703
# round06/round10 scan datasets and are read run-by-run; they have no separate
# sample name in the BESIII dataset table.
# =============================================================================
scan_points = [
  DatasetManager.real_data.find("703_3900"),     # 3.8962 GeV,   52.6 pb^-1
  DatasetManager.real_data.find("703_4009"),     # 4.0076 GeV,  482.0 pb^-1
  DatasetManager.real_data.find("703_4090"),     # 4.0855 GeV,   52.6 pb^-1
  DatasetManager.real_data.find("703_4190scan"), # 4.1886 GeV,   43.1 pb^-1
  DatasetManager.real_data.find("703_4210scan"), # 4.2077 GeV,   54.6 pb^-1
  DatasetManager.real_data.find("703_4220scan"), # 4.2171 GeV,   54.1 pb^-1
  DatasetManager.real_data.find("703_4230"),     # 4.2263 GeV, 1091.7 pb^-1
  DatasetManager.real_data.find("703_4245"),     # 4.2417 GeV,   55.6 pb^-1
  DatasetManager.real_data.find("703_4260"),     # 4.2580 GeV,  825.7 pb^-1
  DatasetManager.real_data.find("703_4310"),     # 4.3079 GeV,   44.9 pb^-1
  DatasetManager.real_data.find("703_4360"),     # 4.3583 GeV,  539.8 pb^-1
  DatasetManager.real_data.find("703_4390"),     # 4.3874 GeV,   55.2 pb^-1
  DatasetManager.real_data.find("703_4420"),     # 4.4156 GeV, 1073.6 pb^-1
  DatasetManager.real_data.find("703_4470"),     # 4.4671 GeV,  109.9 pb^-1
  DatasetManager.real_data.find("703_4530"),     # 4.5271 GeV,  110.0 pb^-1
  DatasetManager.real_data.find("703_4575"),     # 4.5745 GeV,   47.7 pb^-1
  DatasetManager.real_data.find("703_4600")      # 4.5995 GeV,  566.9 pb^-1
]

# Inclusive MC is produced only for a subset of the energy points.
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4009"),
  DatasetManager.inclusive_mc.find("703_4190scan"),
  DatasetManager.inclusive_mc.find("703_4210scan"),
  DatasetManager.inclusive_mc.find("703_4220scan"),
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4600")
]

# =============================================================================
# Mode 01 : eta_c -> p pbar
# =============================================================================
decay_card_m01 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 p+ anti-p-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m01 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_ppbar"
  config.events        = 40_000
  config.decay_card    = decay_card_m01
  config.cross_section = :default
end

alg_m01 = Algorithm.new("PiPiHcEtacPPbar")
alg_m01.set_header(["PiPiHcEtacPPbarAlg/PiPiHcEtacPPbar.h"])
# Multi-energy cross-section measurement: ECMS is injected per job at run time,
# so it is deliberately NOT set via set_constant here.
sel_m01 = Selection.new
sel_m01.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=1"    # p
         nChrn     ">=1"    # pbar
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"   # the E1 photon of h_c -> gamma eta_c
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :proton, against: [:kaon, :pion]
         nprp ">=1"; nprm ">=1"
       }
       .kinematic_fit([:pip, :pim, :gamma, :prp, :prm]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m01.note(:kinematic_fit_chi2, "The analysis applies a 4C kinematic fit under the e+e- -> pi+ pi- gamma X_i hypothesis; the reference analysis (Ref. [17] of the Letter) requires chi2_4C < 25. The loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "If more than one candidate satisfies the selection, the one with the smallest chi2_4C is retained (handled by the minimum-chi2 combination selection of the kinematic fit).")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window; applied on the fit-updated four-momenta in ROOT.")
      .note(:eta_c_mass_window, "The eta_c candidate is selected with a mass window of about +-50 MeV/c^2 around the nominal eta_c mass (efficiency about 84%) for final states with only charged tracks or K_S0; applied in ROOT.")
      .note(:hc_signal_region, "The h_c signal region is [3.515, 3.535] GeV/c^2 for M(gamma eta_c); the two sideband regions [3.475, 3.495] and [3.555, 3.575] GeV/c^2 and the scale factor f = 0.5 are used for the R-scan sample. All of this is applied in ROOT.")
      .note(:pid_selection, "Charged-particle identification uses dE/dx and TOF information; the proton hypothesis is required for this mode.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), with B1 = B(h_c -> gamma eta_c), eps_i the efficiency and B2(i) the branching fraction of the i-th eta_c decay mode; (1+delta) is the ISR correction factor computed in QED and |1+Pi|^2 the vacuum-polarisation factor.")
      .note(:signal_extraction, "For the high-luminosity points the 16 eta_c decay modes are fitted simultaneously to the M(gamma eta_c) spectrum with the mode yields constrained by the corresponding branching fractions; elsewhere the spectrum summed over all modes is fitted; for the R-scan sample the signal is obtained by counting: n_hc^obs = n^sig - f n^side. Performed in ROOT.")
alg_m01.with_decay_card(decay_card_m01).apply(sel_m01)
alg_m01.execute_on(scan_points + incMC_points + exMC_m01)

# =============================================================================
# Mode 02 : eta_c -> 2 (pi+ pi-)
# =============================================================================
decay_card_m02 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m02 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_2pipi"
  config.events        = 40_000
  config.decay_card    = decay_card_m02
  config.cross_section = :default
end

alg_m02 = Algorithm.new("PiPiHcEtac2PiPi")
alg_m02.set_header(["PiPiHcEtac2PiPiAlg/PiPiHcEtac2PiPi.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m02 = Selection.new
sel_m02.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=2"
         nChrn     ">=2"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"   # the E1 photon of h_c -> gamma eta_c
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]
       }
       .kinematic_fit([:pip, :pim, :gamma, :pip, :pim, :pip, :pim]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m02.note(:kinematic_fit_chi2, "The analysis applies a 4C kinematic fit under the e+e- -> pi+ pi- gamma X_i hypothesis; the reference analysis (Ref. [17]) requires chi2_4C < 25. The loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-50 MeV/c^2 around the nominal eta_c mass (efficiency about 84%) for this charged-only final state (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 for the R-scan sample (ROOT level).")
      .note(:pid_selection, "Charged-particle identification uses dE/dx and TOF information; all charged tracks are required to be consistent with the pion hypothesis.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), with B1 = B(h_c -> gamma eta_c); (1+delta) and |1+Pi|^2 obtained from QED calculations, iterated.")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes for the high-luminosity points, summed spectrum elsewhere, and sideband counting for the R-scan sample (ROOT level).")
alg_m02.with_decay_card(decay_card_m02).apply(sel_m02)
alg_m02.execute_on(scan_points + incMC_points + exMC_m02)

# =============================================================================
# Mode 03 : eta_c -> 2 (K+ K-)
# =============================================================================
decay_card_m03 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- K+ K-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m03 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_2KK"
  config.events        = 40_000
  config.decay_card    = decay_card_m03
  config.cross_section = :default
end

alg_m03 = Algorithm.new("PiPiHcEtac2KK")
alg_m03.set_header(["PiPiHcEtac2KKAlg/PiPiHcEtac2KK.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m03 = Selection.new
sel_m03.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=2"
         nChrn     ">=2"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :kaon, against: [:pion, :proton]
         nkp ">=2"; nkm ">=2"
       }
       .kinematic_fit([:pip, :pim, :gamma, :kp, :km, :kp, :km]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m03.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-50 MeV/c^2 around the nominal eta_c mass (about 84% efficiency) for this charged-only final state (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; both K+ and K- are required.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m03.with_decay_card(decay_card_m03).apply(sel_m03)
alg_m03.execute_on(scan_points + incMC_points + exMC_m03)

# =============================================================================
# Mode 04 : eta_c -> K+ K- pi+ pi-
# =============================================================================
decay_card_m04 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m04 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_KKPiPi"
  config.events        = 40_000
  config.decay_card    = decay_card_m04
  config.cross_section = :default
end

alg_m04 = Algorithm.new("PiPiHcEtacKKPiPi")
alg_m04.set_header(["PiPiHcEtacKKPiPiAlg/PiPiHcEtacKKPiPi.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m04 = Selection.new
sel_m04.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=2"
         nChrn     ">=2"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :kaon, against: [:pion, :proton]
         identify :pion, against: [:kaon, :proton]
         nkp ">=1"; nkm ">=1"
         npip ">=1"; npim ">=1"
       }
       .kinematic_fit([:pip, :pim, :gamma, :kp, :km, :pip, :pim]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m04.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-50 MeV/c^2 around the nominal eta_c mass (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF combined with the kinematic fit, which resolves the pi/K separation.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m04.with_decay_card(decay_card_m04).apply(sel_m04)
alg_m04.execute_on(scan_points + incMC_points + exMC_m04)

# =============================================================================
# Mode 05 : eta_c -> p pbar pi+ pi-
# =============================================================================
decay_card_m05 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 p+ anti-p- pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m05 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_ppbarpipi"
  config.events        = 40_000
  config.decay_card    = decay_card_m05
  config.cross_section = :default
end

alg_m05 = Algorithm.new("PiPiHcEtacPPbarPiPi")
alg_m05.set_header(["PiPiHcEtacPPbarPiPiAlg/PiPiHcEtacPPbarPiPi.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m05 = Selection.new
sel_m05.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=2"
         nChrn     ">=2"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :proton, against: [:kaon, :pion]
         identify :pion,   against: [:kaon, :proton]
         nprp ">=1"; nprm ">=1"
         npip ">=1"; npim ">=1"
       }
       .kinematic_fit([:pip, :pim, :gamma, :prp, :prm, :pip, :pim]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m05.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-50 MeV/c^2 around the nominal eta_c mass (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; one p, one pbar, one pi+ and one pi- are required.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m05.with_decay_card(decay_card_m05).apply(sel_m05)
alg_m05.execute_on(scan_points + incMC_points + exMC_m05)

# =============================================================================
# Mode 06 : eta_c -> 3 (pi+ pi-)
# =============================================================================
decay_card_m06 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi+ pi- pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m06 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_3pipi"
  config.events        = 40_000
  config.decay_card    = decay_card_m06
  config.cross_section = :default
end

alg_m06 = Algorithm.new("PiPiHcEtac3PiPi")
alg_m06.set_header(["PiPiHcEtac3PiPiAlg/PiPiHcEtac3PiPi.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m06 = Selection.new
sel_m06.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=3"
         nChrn     ">=3"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]
       }
       .kinematic_fit([:pip, :pim, :gamma, :pip, :pim, :pip, :pim, :pip, :pim]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m06.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-50 MeV/c^2 around the nominal eta_c mass (about 84% efficiency) for this charged-only final state (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; all charged tracks are pions.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m06.with_decay_card(decay_card_m06).apply(sel_m06)
alg_m06.execute_on(scan_points + incMC_points + exMC_m06)

# =============================================================================
# Mode 07 : eta_c -> 2 (pi+ pi-) K+ K-
# =============================================================================
decay_card_m07 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- pi+ pi- pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m07 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_KK2pipi"
  config.events        = 40_000
  config.decay_card    = decay_card_m07
  config.cross_section = :default
end

alg_m07 = Algorithm.new("PiPiHcEtacKK2PiPi")
alg_m07.set_header(["PiPiHcEtacKK2PiPiAlg/PiPiHcEtacKK2PiPi.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m07 = Selection.new
sel_m07.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=3"
         nChrn     ">=3"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :kaon, against: [:pion, :proton]
         identify :pion, against: [:kaon, :proton]
         nkp ">=1"; nkm ">=1"
         npip ">=1"; npim ">=1"
       }
       .kinematic_fit([:pip, :pim, :gamma, :kp, :km, :pip, :pim, :pip, :pim]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m07.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-50 MeV/c^2 around the nominal eta_c mass (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF combined with the kinematic fit, which resolves the pi/K separation.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m07.with_decay_card(decay_card_m07).apply(sel_m07)
alg_m07.execute_on(scan_points + incMC_points + exMC_m07)

# =============================================================================
# Mode 08 : eta_c -> K_S0 K+- pi-+  (charge-conjugate mode)
# =============================================================================
decay_card_m08 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K_S0 K+ pi-   PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m08 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_KSKpi"
  config.events        = 40_000
  config.decay_card    = decay_card_m08
  config.cross_section = :default
end

alg_m08 = Algorithm.new("PiPiHcEtacKSKPi")
alg_m08.set_header(["PiPiHcEtacKSKPiAlg/PiPiHcEtacKSKPi.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m08 = Selection.new
sel_m08.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=2"    # K+ and the pi+ of K_S0
         nChrn     ">=2"    # pi- and the pi- of K_S0
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :kaon, against: [:pion, :proton]
         nkp ">=1"; nkm ">=1"
       }
       .secondary_vertex_fit([:pip, :pim]) {
         build_virtual_particle(:K_S0).by_minimizing_mass_difference
         remove_used_particle_from_candidate_list
       }
       .kinematic_fit([:pip, :pim, :gamma, :K_S0, :kp, :pim]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m08.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:charge_conjugate_mode, "The K_S0 K+- pi-+ final state contains both charge-conjugate combinations; the decay card and the reconstruction are written for K_S0 K+ pi-, the charge-conjugate assignment being covered by the charge-symmetric selection.")
      .note(:ks0_selection, "K_S0 -> pi+ pi- candidates are reconstructed with a secondary-vertex fit constraining the two tracks to a common decay vertex, with the invariant mass required to be consistent with the nominal K_S0 mass (about +-20 MeV/c^2). The 1.2% per-K_S0 reconstruction systematic is quoted in the Letter.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-50 MeV/c^2 around the nominal eta_c mass (about 84% efficiency) for K_S0 final states (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; the K+ is separated from the pions with the kinematic fit.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m08.with_decay_card(decay_card_m08).apply(sel_m08)
alg_m08.execute_on(scan_points + incMC_points + exMC_m08)

# =============================================================================
# Mode 09 : eta_c -> K_S0 K+- pi-+ pi+ pi-
# =============================================================================
decay_card_m09 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K_S0 K+ pi- pi+ pi-   PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m09 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_KSK2pi"
  config.events        = 40_000
  config.decay_card    = decay_card_m09
  config.cross_section = :default
end

alg_m09 = Algorithm.new("PiPiHcEtacKSK2Pi")
alg_m09.set_header(["PiPiHcEtacKSK2PiAlg/PiPiHcEtacKSK2Pi.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m09 = Selection.new
sel_m09.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=3"    # K+, pi+ of K_S0 and pi+
         nChrn     ">=3"    # pi-, pi- of K_S0 and pi-
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :kaon, against: [:pion, :proton]
         nkp ">=1"; nkm ">=1"
       }
       .secondary_vertex_fit([:pip, :pim]) {
         build_virtual_particle(:K_S0).by_minimizing_mass_difference
         remove_used_particle_from_candidate_list
       }
       .kinematic_fit([:pip, :pim, :gamma, :K_S0, :kp, :pim, :pip, :pim]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m09.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:charge_conjugate_mode, "The K_S0 K+- pi-+ pi+ pi- final state contains both charge-conjugate combinations; the reconstruction is written for K_S0 K+ pi- pi+ pi-.")
      .note(:ks0_selection, "K_S0 -> pi+ pi- candidates are reconstructed with a common decay-vertex fit and a mass requirement consistent with the nominal K_S0 mass (about +-20 MeV/c^2).")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-50 MeV/c^2 around the nominal eta_c mass (about 84% efficiency) for K_S0 final states (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; the K+ is separated from the pions with the kinematic fit.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m09.with_decay_card(decay_card_m09).apply(sel_m09)
alg_m09.execute_on(scan_points + incMC_points + exMC_m09)

# =============================================================================
# Mode 10 : eta_c -> K+ K- pi0
# =============================================================================
decay_card_m10 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- pi0   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m10 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_KKpi0"
  config.events        = 40_000
  config.decay_card    = decay_card_m10
  config.cross_section = :default
end

alg_m10 = Algorithm.new("PiPiHcEtacKKPi0")
alg_m10.set_header(["PiPiHcEtacKKPi0Alg/PiPiHcEtacKKPi0.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m10 = Selection.new
sel_m10.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=1"    # K+
         nChrn     ">=1"    # K-
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=3"   # E1 photon + the two photons of the pi0
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :kaon, against: [:pion, :proton]
         nkp ">=1"; nkm ">=1"
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
         chi2_cut 25
         npi0 ">=1"
       }
       .kinematic_fit([:pip, :pim, :gamma, :kp, :km, :pi0]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m10.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:pi0_mass_window, "pi0 candidates are formed from photon pairs with an invariant mass consistent with the nominal pi0 mass (about +-15 MeV/c^2), improved by the 1C mass-constrained fit (kalman_kinematic_fit); the 1% per-pi0 systematic from the mass window is quoted in the Letter.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-45 MeV/c^2 around the nominal eta_c mass (about 80% efficiency) for final states containing a pi0 or eta (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; kaon hypotheses are required, the pi0 being identified through its gamma gamma decay.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m10.with_decay_card(decay_card_m10).apply(sel_m10)
alg_m10.execute_on(scan_points + incMC_points + exMC_m10)

# =============================================================================
# Mode 11 : eta_c -> p pbar pi0
# =============================================================================
decay_card_m11 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 p+ anti-p- pi0   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m11 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_ppbarpi0"
  config.events        = 40_000
  config.decay_card    = decay_card_m11
  config.cross_section = :default
end

alg_m11 = Algorithm.new("PiPiHcEtacPPbarPi0")
alg_m11.set_header(["PiPiHcEtacPPbarPi0Alg/PiPiHcEtacPPbarPi0.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m11 = Selection.new
sel_m11.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=1"
         nChrn     ">=1"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=3"   # E1 photon + the two photons of the pi0
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :proton, against: [:kaon, :pion]
         nprp ">=1"; nprm ">=1"
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
         chi2_cut 25
         npi0 ">=1"
       }
       .kinematic_fit([:pip, :pim, :gamma, :prp, :prm, :pi0]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m11.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:pi0_mass_window, "pi0 candidates: photon pairs with an invariant mass consistent with the nominal pi0 mass, improved by a 1C mass-constrained fit.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-45 MeV/c^2 around the nominal eta_c mass (about 80% efficiency) for this pi0 mode (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; the proton and antiproton hypotheses are required.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m11.with_decay_card(decay_card_m11).apply(sel_m11)
alg_m11.execute_on(scan_points + incMC_points + exMC_m11)

# =============================================================================
# Mode 12 : eta_c -> K+ K- eta
# =============================================================================
decay_card_m12 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- eta   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m12 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_KKeta"
  config.events        = 40_000
  config.decay_card    = decay_card_m12
  config.cross_section = :default
end

alg_m12 = Algorithm.new("PiPiHcEtacKKEta")
alg_m12.set_header(["PiPiHcEtacKKEtaAlg/PiPiHcEtacKKEta.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m12 = Selection.new
sel_m12.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=1"
         nChrn     ">=1"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=3"   # E1 photon + the two photons of the eta
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :kaon, against: [:pion, :proton]
         nkp ">=1"; nkm ">=1"
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
         chi2_cut 25
         neta ">=1"
       }
       .kinematic_fit([:pip, :pim, :gamma, :kp, :km, :eta]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m12.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:eta_mass_window, "eta candidates: photon pairs with an invariant mass consistent with the nominal eta mass, improved by a 1C mass-constrained fit; the 1% per-eta systematic from the mass window is quoted in the Letter.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-45 MeV/c^2 around the nominal eta_c mass (about 80% efficiency) for final states containing a pi0 or eta (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; kaon hypotheses are required, the eta being identified through its gamma gamma decay.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m12.with_decay_card(decay_card_m12).apply(sel_m12)
alg_m12.execute_on(scan_points + incMC_points + exMC_m12)

# =============================================================================
# Mode 13 : eta_c -> pi+ pi- eta
# =============================================================================
decay_card_m13 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- eta   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m13 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_pipieta"
  config.events        = 40_000
  config.decay_card    = decay_card_m13
  config.cross_section = :default
end

alg_m13 = Algorithm.new("PiPiHcEtacPiPiEta")
alg_m13.set_header(["PiPiHcEtacPiPiEtaAlg/PiPiHcEtacPiPiEta.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m13 = Selection.new
sel_m13.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=1"
         nChrn     ">=1"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=3"   # E1 photon + the two photons of the eta
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
         chi2_cut 25
         neta ">=1"
       }
       .kinematic_fit([:pip, :pim, :gamma, :pip, :pim, :eta]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m13.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:eta_mass_window, "eta candidates: photon pairs with an invariant mass consistent with the nominal eta mass, improved by a 1C mass-constrained fit.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-45 MeV/c^2 around the nominal eta_c mass (about 80% efficiency) for this eta mode (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; the charged tracks are pions, the eta being identified through its gamma gamma decay.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m13.with_decay_card(decay_card_m13).apply(sel_m13)
alg_m13.execute_on(scan_points + incMC_points + exMC_m13)

# =============================================================================
# Mode 14 : eta_c -> 2 (pi+ pi-) eta
# =============================================================================
decay_card_m14 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi+ pi- eta   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m14 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_2pipieta"
  config.events        = 40_000
  config.decay_card    = decay_card_m14
  config.cross_section = :default
end

alg_m14 = Algorithm.new("PiPiHcEtac2PiPiEta")
alg_m14.set_header(["PiPiHcEtac2PiPiEtaAlg/PiPiHcEtac2PiPiEta.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m14 = Selection.new
sel_m14.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=2"
         nChrn     ">=2"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=3"   # E1 photon + the two photons of the eta
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
         chi2_cut 25
         neta ">=1"
       }
       .kinematic_fit([:pip, :pim, :gamma, :pip, :pim, :pip, :pim, :eta]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m14.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained.")
      .note(:eta_mass_window, "eta candidates: photon pairs with an invariant mass consistent with the nominal eta mass, improved by a 1C mass-constrained fit.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-45 MeV/c^2 around the nominal eta_c mass (about 80% efficiency) for this eta mode (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; all charged tracks are pions, the eta being identified through its gamma gamma decay.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m14.with_decay_card(decay_card_m14).apply(sel_m14)
alg_m14.execute_on(scan_points + incMC_points + exMC_m14)

# =============================================================================
# Mode 15 : eta_c -> pi+ pi- pi0 pi0
# =============================================================================
decay_card_m15 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi0 pi0   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m15 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_pipi2pi0"
  config.events        = 40_000
  config.decay_card    = decay_card_m15
  config.cross_section = :default
end

alg_m15 = Algorithm.new("PiPiHcEtacPiPi2Pi0")
alg_m15.set_header(["PiPiHcEtacPiPi2Pi0Alg/PiPiHcEtacPiPi2Pi0.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m15 = Selection.new
sel_m15.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=1"
         nChrn     ">=1"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=5"   # E1 photon + the four photons of the two pi0
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
         chi2_cut 25
         npi0 ">=2"
       }
       .kinematic_fit([:pip, :pim, :gamma, :pip, :pim, :pi0, :pi0]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m15.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained; the pi0 pairing ambiguity is resolved together with the fit.")
      .note(:pi0_mass_window, "pi0 candidates are formed from photon pairs with an invariant mass consistent with the nominal pi0 mass, improved by the 1C mass-constrained fit; two pi0 are required.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-45 MeV/c^2 around the nominal eta_c mass (about 80% efficiency) for this pi0 mode (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; the charged tracks are pions.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m15.with_decay_card(decay_card_m15).apply(sel_m15)
alg_m15.execute_on(scan_points + incMC_points + exMC_m15)

# =============================================================================
# Mode 16 : eta_c -> 2 (pi+ pi- pi0)
# =============================================================================
decay_card_m16 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi0 pi+ pi- pi0   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m16 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_2pipipi0"
  config.events        = 40_000
  config.decay_card    = decay_card_m16
  config.cross_section = :default
end

alg_m16 = Algorithm.new("PiPiHcEtac2PiPiPi0")
alg_m16.set_header(["PiPiHcEtac2PiPiPi0Alg/PiPiHcEtac2PiPiPi0.h"])
# Multi-energy cross-section measurement: ECMS is deliberately NOT set.
sel_m16 = Selection.new
sel_m16.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     ">=2"
         nChrn     ">=2"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=5"   # E1 photon + the four photons of the two pi0
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
         chi2_cut 25
         npi0 ">=2"
       }
       .kinematic_fit([:pip, :pim, :gamma, :pip, :pim, :pi0, :pip, :pim, :pi0]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_m16.note(:kinematic_fit_chi2, "4C kinematic fit under e+e- -> pi+ pi- gamma X_i; the reference analysis requires chi2_4C < 25, applied in ROOT, the loose default (200) being used in BOSS.")
      .note(:combination_selection, "The candidate with the smallest chi2_4C is retained; the pi0 pairing ambiguity is resolved together with the fit.")
      .note(:pi0_mass_window, "pi0 candidates are formed from photon pairs with an invariant mass consistent with the nominal pi0 mass, improved by the 1C mass-constrained fit; two pi0 are required.")
      .note(:e1_photon_energy, "The photon from h_c -> gamma eta_c is required to lie in an optimised energy window (ROOT level).")
      .note(:eta_c_mass_window, "Mass window of about +-45 MeV/c^2 around the nominal eta_c mass (about 80% efficiency) for this pi0 mode (ROOT level).")
      .note(:hc_signal_region, "h_c signal region [3.515, 3.535] GeV/c^2 in M(gamma eta_c); sidebands [3.475, 3.495] and [3.555, 3.575] GeV/c^2 with f = 0.5 (ROOT level).")
      .note(:pid_selection, "Charged-particle identification from dE/dx and TOF; the charged tracks are pions.")
      .note(:born_cross_section, "sigma^B = n_hc^obs / (L (1+delta) |1+Pi|^2 B1 sum_i eps_i B2(i)), B1 = B(h_c -> gamma eta_c).")
      .note(:signal_extraction, "Simultaneous fit of the 16 eta_c modes (high-luminosity points), summed spectrum elsewhere, sideband counting for the R-scan sample (ROOT level).")
alg_m16.with_decay_card(decay_card_m16).apply(sel_m16)
alg_m16.execute_on(scan_points + incMC_points + exMC_m16)
