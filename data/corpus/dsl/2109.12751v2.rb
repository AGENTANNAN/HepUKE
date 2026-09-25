# Paper: 2109.12751v2
# Title: Cross sections for e+e- -> K+K-pi+pi-(pi0), K+K-K+K-(pi0), pi+pi-pi+pi-(pi0), p pbar pi+pi-(pi0)
#         at sqrt(s) = 3.773-4.600 GeV
#
# 8 independent channels (Rule T1 — each has its own Algorithm):
#   Ch1: e+e- -> K+ K- pi+ pi-        (4C fit, chi2<50)
#   Ch2: e+e- -> K+ K- K+ K-          (4C fit, chi2<50)
#   Ch3: e+e- -> pi+ pi- pi+ pi-      (4C fit, chi2<50)
#   Ch4: e+e- -> p pbar pi+ pi-       (4C fit, chi2<50)
#   Ch5: e+e- -> K+ K- pi+ pi- pi0    (5C fit, chi2<50)
#   Ch6: e+e- -> K+ K- K+ K- pi0      (5C fit, chi2<50)
#   Ch7: e+e- -> pi+ pi- pi+ pi- pi0  (5C fit, chi2<50)
#   Ch8: e+e- -> p pbar pi+ pi- pi0   (5C fit, chi2<50)
#
# ConExc generator for continuum production (dressed cross section vs sqrt(s)).
# Data: 40 energy points from 3.773 to 4.600 GeV (BOSS 703/706/707).

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ============================================================
# Energy points — representative set covering the full range
# BOSS 703 points: 4.009, 4.178-4.600 GeV region
# BOSS 705 points: 4.130, 4.160 GeV
# BOSS 706 points: 4.612, 4.628, 4.640, 4.660, 4.680, 4.700 GeV
# BOSS 707 points: 4.740, 4.750, 4.780, 4.840, 4.914 GeV
# plus 3.773 (psi(3770)) and 3.808, 3.867, 3.871, 3.896 GeV (BOSS 703)
# ============================================================

# Representative energy points (full list would include all 40):
scan_points = [
  DatasetManager.real_data.find("712_3773"),
  DatasetManager.real_data.find("703_4009"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4600"),
]

scan_incMC = scan_points.map { |dp| DatasetManager.inclusive_mc.find(dp.sample_name) rescue nil }.compact

# ConExc decay card for e+e- -> K+ K- pi+ pi- (mode 14)
# Omit 'Particle vpho' — DSL auto-injects per energy point
conexc_kkpipi = <<~DECAYCARD
  Decay vpho
  1 ConExc 14;
  Enddecay
  Decay vhdr
  1 K+ K- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

conexc_kkkk = <<~DECAYCARD
  Decay vpho
  1 ConExc 16;
  Enddecay
  Decay vhdr
  1 K+ K- K+ K- PHSP;
  Enddecay
  End
DECAYCARD

conexc_4pi = <<~DECAYCARD
  Decay vpho
  1 ConExc 12;
  Enddecay
  Decay vhdr
  1 pi+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

conexc_ppbar_pipi = <<~DECAYCARD
  Decay vpho
  1 ConExc -2 2212 -2112 211 -211;   # p pbar pi+ pi- (form D: user xs)
  Enddecay
  End
DECAYCARD

conexc_kkpipi_pi0 = <<~DECAYCARD
  Decay vpho
  1 ConExc 19;
  Enddecay
  Decay vhdr
  1 K+ K- pi+ pi- pi0 PHSP;
  Enddecay
  End
DECAYCARD

conexc_kkkk_pi0 = <<~DECAYCARD
  Decay vpho
  1 ConExc -2 321 -321 321 -321 111;   # K+ K- K+ K- pi0 (form D)
  Enddecay
  End
DECAYCARD

conexc_4pi_pi0 = <<~DECAYCARD
  Decay vpho
  1 ConExc 17;
  Enddecay
  Decay vhdr
  1 pi+ pi- pi+ pi- pi0 PHSP;
  Enddecay
  End
DECAYCARD

conexc_ppbar_pipi_pi0 = <<~DECAYCARD
  Decay vpho
  1 ConExc -2 2212 -2112 211 -211 111;   # p pbar pi+ pi- pi0 (form D)
  Enddecay
  End
DECAYCARD

# Signal MC — use create_exclusive_mc_for for scan points
sig_kkpipi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_conexc_kkpipi"
  config.events        = 100_000
  config.decay_card    = conexc_kkpipi
  config.cross_section = :default
end

sig_kkkk = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_conexc_kkkk"
  config.events        = 100_000
  config.decay_card    = conexc_kkkk
  config.cross_section = :default
end

sig_4pi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_conexc_4pi"
  config.events        = 100_000
  config.decay_card    = conexc_4pi
  config.cross_section = :default
end

sig_ppbar_pipi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_conexc_ppbarpipi"
  config.events        = 100_000
  config.decay_card    = conexc_ppbar_pipi
  config.cross_section = :default
end

sig_kkpipi_pi0 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_conexc_kkpipipi0"
  config.events        = 100_000
  config.decay_card    = conexc_kkpipi_pi0
  config.cross_section = :default
end

sig_kkkk_pi0 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_conexc_kkkkpi0"
  config.events        = 100_000
  config.decay_card    = conexc_kkkk_pi0
  config.cross_section = :default
end

sig_4pi_pi0 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_conexc_4pipi0"
  config.events        = 100_000
  config.decay_card    = conexc_4pi_pi0
  config.cross_section = :default
end

sig_ppbar_pipi_pi0 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_conexc_ppbarpipipi0"
  config.events        = 100_000
  config.decay_card    = conexc_ppbar_pipi_pi0
  config.cross_section = :default
end

# ============================================================
# Algorithm 1: e+e- -> K+ K- pi+ pi-  (4C fit)
# Charged tracks: 4 (K+, K-, pi+, pi-)
# PID: K/pi separation
# Vetoes: |M_pipi - M_Jpsi|>15 MeV, |M_Kpi - M_D|>20 MeV
# E_EMC/p < 0.8, cos_theta_pipi < 0.9
# ============================================================

alg_kkpipi = Algorithm.new("EEtoKKPiPi")
alg_kkpipi.set_header(["EEtoKKPiPiAlg/EEtoKKPiPi.h"])
          .set_constant({ "ECMS" => [:double, 4.226] })

sel_kkpipi = Selection.new
  .select_track { nChrp "==2"; nChrn "==2"; cos_theta 0.93; Vz 10.0; Vr 1.0 }
  .select_photon { nGam ">=0" }
  .pid(method: :probability) do
    identify :kp, against: :pion
    identify :km, against: :pion
    identify :pip, against: :kaon
    identify :pim, against: :kaon
    prob_cut 0.001
  end
  .kinematic_fit([:kp, :km, :pip, :pim]) do
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

alg_kkpipi.note(:background_veto, "|M(pipi) - M_Jpsi| > 15 MeV/c2. |M(Kpi) - M_D0| > 20 MeV/c2. " \
                                   "E_EMC/p < 0.8 for each track. cos_theta(pipi) < 0.9.")
alg_kkpipi.note(:cross_section, "Dressed cross section: sigma = (N_obs - N_bkg) / (L * epsilon0 * kappa). " \
                                 "ISR correction factor kappa from ConExc. " \
                                 "Non-peaking background from chi2 fit. Peaking background < 0.8%.")
alg_kkpipi.with_decay_card(conexc_kkpipi).apply(sel_kkpipi)
alg_kkpipi.execute_on(scan_points + scan_incMC + sig_kkpipi)

# ============================================================
# Algorithm 2: e+e- -> K+ K- K+ K-  (4C fit)
# Charged tracks: 4 (K+, K-, K+, K-)
# PID: K/pi separation
# ============================================================

alg_kkkk = Algorithm.new("EEtoKKKK")
alg_kkkk.set_header(["EEtoKKKKAlg/EEtoKKKK.h"])
        .set_constant({ "ECMS" => [:double, 4.226] })

sel_kkkk = Selection.new
  .select_track { nChrp "==2"; nChrn "==2"; cos_theta 0.93; Vz 10.0; Vr 1.0 }
  .select_photon { nGam ">=0" }
  .pid(method: :probability) do
    identify :kp, against: :pion
    identify :km, against: :pion
    prob_cut 0.001
  end
  .kinematic_fit([:kp, :km, :kp, :km]) do
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

alg_kkkk.note(:background_veto, "E_EMC/p < 0.8 for each track. " \
                                 "cos_theta(pipi) checks for gamma conversion.")
alg_kkkk.note(:cross_section, "Dressed cross section measured at each energy point. " \
                               "Continuum + possible resonance contributions.")
alg_kkkk.with_decay_card(conexc_kkkk).apply(sel_kkkk)
alg_kkkk.execute_on(scan_points + scan_incMC + sig_kkkk)

# ============================================================
# Algorithm 3: e+e- -> pi+ pi- pi+ pi-  (4C fit)
# Charged tracks: 4 (pi+, pi-, pi+, pi-)
# Vetoes: |M_pipi - M_Jpsi|>15 MeV, |M_4pi - M_psi(2S)|>20 MeV
# ============================================================

alg_4pi = Algorithm.new("EEto4Pi")
alg_4pi.set_header(["EEto4PiAlg/EEto4Pi.h"])
       .set_constant({ "ECMS" => [:double, 4.226] })

sel_4pi = Selection.new
  .select_track { nChrp "==2"; nChrn "==2"; cos_theta 0.93; Vz 10.0; Vr 1.0 }
  .select_photon { nGam ">=0" }
  .pid(method: :probability) do
    identify :pip, against: :kaon
    identify :pim, against: :kaon
    prob_cut 0.001
  end
  .kinematic_fit([:pip, :pim, :pip, :pim]) do
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

alg_4pi.note(:background_veto, "|M(pipi) - M_Jpsi| > 15 MeV/c2. |M(4pi) - M_psi(2S)| > 20 MeV/c2. " \
                                "E_EMC/p < 0.8. cos_theta(pipi) < 0.9.")
alg_4pi.note(:psi4040, "Evidence for psi(4040) -> pi+ pi- pi+ pi- pi0 with 3.6 sigma significance. " \
                        "This is in the 5pi channel (Algorithm 7), not this one.")
alg_4pi.with_decay_card(conexc_4pi).apply(sel_4pi)
alg_4pi.execute_on(scan_points + scan_incMC + sig_4pi)

# ============================================================
# Algorithm 4: e+e- -> p pbar pi+ pi-  (4C fit)
# Charged tracks: 4 (p, pbar, pi+, pi-)
# PID: proton/pion separation
# Vetoes: |M_ppbar - M_Jpsi|>15 MeV, |M_ppbar - M_psi(2S)|>20 MeV
# ============================================================

alg_ppbar_pipi = Algorithm.new("EEtoPPbarPiPi")
alg_ppbar_pipi.set_header(["EEtoPPbarPiPiAlg/EEtoPPbarPiPi.h"])
              .set_constant({ "ECMS" => [:double, 4.226] })

sel_ppbar_pipi = Selection.new
  .select_track { nChrp "==2"; nChrn "==2"; cos_theta 0.93; Vz 10.0; Vr 1.0 }
  .select_photon { nGam ">=0" }
  .pid(method: :probability) do
    identify :prp, against: [:pion, :kaon]
    identify :prm, against: [:pion, :kaon]
    identify :pip, against: :proton
    identify :pim, against: :proton
    prob_cut 0.001
  end
  .kinematic_fit([:prp, :prm, :pip, :pim]) do
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

alg_ppbar_pipi.note(:background_veto, "|M(ppbar) - M_Jpsi| > 15 MeV/c2. |M(ppbar) - M_psi(2S)| > 20 MeV/c2. " \
                                       "E_EMC/p < 0.8 for each track. cos_theta(pipi) < 0.9.")
alg_ppbar_pipi.note(:conexc_form_d, "p pbar pi+ pi- uses ConExc form D (-2) with user cross-section table xs_user.txt.")
alg_ppbar_pipi.with_decay_card(conexc_ppbar_pipi).apply(sel_ppbar_pipi)
alg_ppbar_pipi.execute_on(scan_points + scan_incMC + sig_ppbar_pipi)

# ============================================================
# Algorithm 5: e+e- -> K+ K- pi+ pi- pi0  (5C fit)
# Charged tracks: 4, photons: 2 (from pi0)
# PID: K/pi separation
# 5C = 4C + pi0 mass constraint
# Vetoes: |M_KK - M_chic0|>50 MeV, |M_Kpi - M_D|>20 MeV, |M_pipi - M_Ks|>30 MeV
# ============================================================

alg_kkpipi_pi0 = Algorithm.new("EEtoKKPiPiPi0")
alg_kkpipi_pi0.set_header(["EEtoKKPiPiPi0Alg/EEtoKKPiPiPi0.h"])
              .set_constant({ "ECMS" => [:double, 4.226] })

sel_kkpipi_pi0 = Selection.new
  .select_track { nChrp "==2"; nChrn "==2"; cos_theta 0.93; Vz 10.0; Vr 1.0 }
  .select_photon { nGam ">=2"; energyThreshold_b 0.025; energyThreshold_e 0.050 }
  .pid(method: :probability) do
    identify :kp, against: :pion
    identify :km, against: :pion
    identify :pip, against: :kaon
    identify :pim, against: :kaon
    prob_cut 0.001
  end
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

alg_kkpipi_pi0.note(:background_veto, "|M(KK) - M_chic0| > 50 MeV/c2. |M(Kpi) - M_D| > 20 MeV/c2. " \
                                       "|M(pipi) - M_Ks0| > 30 MeV/c2. E_EMC/p < 0.8.")
alg_kkpipi_pi0.note(:photon_timing, "EMC time within [0, 700] ns of event start time.")
alg_kkpipi_pi0.with_decay_card(conexc_kkpipi_pi0).apply(sel_kkpipi_pi0)
alg_kkpipi_pi0.execute_on(scan_points + scan_incMC + sig_kkpipi_pi0)

# ============================================================
# Algorithm 6: e+e- -> K+ K- K+ K- pi0  (5C fit)
# Charged tracks: 4 (K+, K-, K+, K-), photons: 2 (from pi0)
# ============================================================

alg_kkkk_pi0 = Algorithm.new("EEtoKKKKPi0")
alg_kkkk_pi0.set_header(["EEtoKKKKPi0Alg/EEtoKKKKPi0.h"])
            .set_constant({ "ECMS" => [:double, 4.226] })

sel_kkkk_pi0 = Selection.new
  .select_track { nChrp "==2"; nChrn "==2"; cos_theta 0.93; Vz 10.0; Vr 1.0 }
  .select_photon { nGam ">=2"; energyThreshold_b 0.025; energyThreshold_e 0.050 }
  .pid(method: :probability) do
    identify :kp, against: :pion
    identify :km, against: :pion
    prob_cut 0.001
  end
  .kinematic_fit([:kp, :km, :kp, :km, :gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

alg_kkkk_pi0.note(:background_veto, "E_EMC/p cuts and gamma conversion vetoes.")
alg_kkkk_pi0.note(:first_measurement, "First measurement of e+e- -> K+K-K+K-pi0 cross sections.")
alg_kkkk_pi0.with_decay_card(conexc_kkkk_pi0).apply(sel_kkkk_pi0)
alg_kkkk_pi0.execute_on(scan_points + scan_incMC + sig_kkkk_pi0)

# ============================================================
# Algorithm 7: e+e- -> pi+ pi- pi+ pi- pi0  (5C fit)
# Charged tracks: 4, photons: 2 (from pi0)
# Vetoes: |M_pipi - M_Jpsi|>15 MeV, |M_pipi - M_chic0|>50 MeV,
#         |M_pipi - M_Ks0|>30 MeV
# Evidence for psi(4040) -> 5pi (3.6 sigma)
# ============================================================

alg_4pi_pi0 = Algorithm.new("EEto4PiPi0")
alg_4pi_pi0.set_header(["EEto4PiPi0Alg/EEto4PiPi0.h"])
           .set_constant({ "ECMS" => [:double, 4.226] })

sel_4pi_pi0 = Selection.new
  .select_track { nChrp "==2"; nChrn "==2"; cos_theta 0.93; Vz 10.0; Vr 1.0 }
  .select_photon { nGam ">=2"; energyThreshold_b 0.025; energyThreshold_e 0.050 }
  .pid(method: :probability) do
    identify :pip, against: :kaon
    identify :pim, against: :kaon
    prob_cut 0.001
  end
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

alg_4pi_pi0.note(:background_veto, "|M(pipi) - M_Jpsi| > 15 MeV/c2. |M(pipi) - M_chic0| > 50 MeV/c2. " \
                                    "|M(pipi) - M_Ks0| > 30 MeV/c2. E_EMC/p < 0.8.")
alg_4pi_pi0.note(:psi4040, "Evidence for psi(4040) -> pi+pi-pi+pi-pi0 with 3.6 sigma. " \
                            "BW fit with M=4039 MeV, Gamma=80 MeV. BF ~ (3.5e-5 or 2.4e-2)% (two solutions).")
alg_4pi_pi0.with_decay_card(conexc_4pi_pi0).apply(sel_4pi_pi0)
alg_4pi_pi0.execute_on(scan_points + scan_incMC + sig_4pi_pi0)

# ============================================================
# Algorithm 8: e+e- -> p pbar pi+ pi- pi0  (5C fit)
# Charged tracks: 4 (p, pbar, pi+, pi-), photons: 2 (from pi0)
# First measurement
# ============================================================

alg_ppbar_pipi_pi0 = Algorithm.new("EEtoPPbarPiPiPi0")
alg_ppbar_pipi_pi0.set_header(["EEtoPPbarPiPiPi0Alg/EEtoPPbarPiPiPi0.h"])
                  .set_constant({ "ECMS" => [:double, 4.226] })

sel_ppbar_pipi_pi0 = Selection.new
  .select_track { nChrp "==2"; nChrn "==2"; cos_theta 0.93; Vz 10.0; Vr 1.0 }
  .select_photon { nGam ">=2"; energyThreshold_b 0.025; energyThreshold_e 0.050 }
  .pid(method: :probability) do
    identify :prp, against: [:pion, :kaon]
    identify :prm, against: [:pion, :kaon]
    identify :pip, against: :proton
    identify :pim, against: :proton
    prob_cut 0.001
  end
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

alg_ppbar_pipi_pi0.note(:background_veto, "E_EMC/p < 0.8. Gamma conversion vetoes.")
alg_ppbar_pipi_pi0.note(:first_measurement, "First measurement of e+e- -> p pbar pi+ pi- pi0 cross sections.")
alg_ppbar_pipi_pi0.with_decay_card(conexc_ppbar_pipi_pi0).apply(sel_ppbar_pipi_pi0)
alg_ppbar_pipi_pi0.execute_on(scan_points + scan_incMC + sig_ppbar_pipi_pi0)

# ============================================================
# Common notes for all 8 channels
# ============================================================
# - Exact 4 charged tracks required.
# - For pi0 channels: E_gamma > 25(50) MeV barrel(endcap).
# - |cos_theta| < 0.80(0.92) barrel(endcap) for photons.
# - EMC time within [0,700] ns.
# - chi2 < 50 for kinematic fit; best combination kept.
# - E_EMC/p < 0.8 for all tracks (suppress e+e- and two-photon bg).
# - cos_theta(pipi) < 0.9 (suppress gamma conversion).
# - Peaking backgrounds from D/Dbar, J/psi decays (MC estimated).
# - Non-peaking bg from chi2 fit, number floated.
# - Amplitude analysis at 7 energy points for efficiency reweighting.
# - ISR correction factor kappa from ConExc+MC convolution.
# - Full 40 energy points span 3.773-4.600 GeV.