# =============================================================================
# BESIII : Observation of a charged charmoniumlike structure Zc(4020) and
#          search for the Zc(3900) in e+e- -> pi+ pi- h_c
# arXiv:1309.1896v2
#
# h_c is reconstructed via h_c -> gamma eta_c, with eta_c -> X_i in 16
# exclusive hadronic final states.  The Born cross sections are measured at
# 13 CM energies between 3.900 and 4.420 GeV.
# BOSS part only : dataset preparation + event selection up to the 4C fit.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets : the 13 CM energy points (BOSS 703, XYZ scan samples)
### ---------------------------------------------------------------------------
scan_points = [
  DatasetManager.real_data.find("703_3900"),   # 3.900 GeV,  52.8  pb^-1
  DatasetManager.real_data.find("703_4009"),   # 4.009 GeV, 482.0  pb^-1
  DatasetManager.real_data.find("703_4090"),   # 4.090 GeV,  51.0  pb^-1
  DatasetManager.real_data.find("703_4190"),   # 4.190 GeV,  43.0  pb^-1
  DatasetManager.real_data.find("703_4210"),   # 4.210 GeV,  54.7  pb^-1
  DatasetManager.real_data.find("703_4220"),   # 4.220 GeV,  54.6  pb^-1
  DatasetManager.real_data.find("703_4230"),   # 4.230 GeV, 1090.0 pb^-1
  DatasetManager.real_data.find("703_4245"),   # 4.245 GeV,  56.0  pb^-1
  DatasetManager.real_data.find("703_4260"),   # 4.260 GeV, 826.8  pb^-1
  DatasetManager.real_data.find("703_4310"),   # 4.310 GeV,  44.9  pb^-1
  DatasetManager.real_data.find("703_4360"),   # 4.360 GeV, 544.5  pb^-1
  DatasetManager.real_data.find("703_4390"),   # 4.390 GeV,  55.1  pb^-1
  DatasetManager.real_data.find("703_4420")    # 4.420 GeV,  44.7  pb^-1
]
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4009"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4360")
]   # inclusive MC available for the main scan points

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
  config.events        = 50_000
  config.decay_card    = decay_card_m01
  config.cross_section = :default
end

alg_m01 = Algorithm.new("HcEtacPPbar")
alg_m01.set_header(["HcEtacPPbarAlg/HcEtacPPbar.h"])
sel_m01 = Selection.new
sel_m01.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=2"    # pi+ from the primary pair and p
             nChrn     ">=2"    # pi- from the primary pair and pbar
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=1"   # gamma from h_c -> gamma eta_c
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :proton, against: [:kaon, :pion]
             identify :pion,   against: [:kaon, :proton]
             nprp ">=1"; nprm ">=1"
             npip ">=1"; npim ">=1"
           }
          .kinematic_fit([:pip, :pim, :gamma, :prp, :prm]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_m01.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 35 for final states with only charged or K_S0 particles (about 80% efficiency) and chi2_4C < 20 for those with pi0 or eta (about 70%). The loose default (200) is used in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Among all combinations the one minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C is chosen (chi2_PID from dE/dx and TOF, chi2_1C from the mass constraint of the two photons of each pi0/eta). In BOSS the combination is resolved by the minimum chi2_4C of the kinematic fit.")
      .note(:recoil_mass_requirements, "At least one combination with M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 is required.")
      .note(:eta_c_mass_window, "The eta_c candidate window is +/-50 MeV/c^2 around the nominal eta_c mass for final states with only charged or K_S0 particles (efficiency about 85%), and +/-45 MeV/c^2 for those with pi0 or eta (about 80%).")
      .note(:ks0_selection, "K_S0 -> pi+ pi- candidates are selected as described in Ref. [13] (secondary vertex).")
      .note(:track_parameter_correction, "A track parameter correction (Ref. [16]) is applied to the MC charged tracks to improve the data/MC agreement of the chi2_4C distribution; the systematic error of the chi2_4C requirement is taken as half of the correction in efficiency.")
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
  config.events        = 50_000
  config.decay_card    = decay_card_m02
  config.cross_section = :default
end

alg_m02 = Algorithm.new("HcEtac2PiPi")
alg_m02.set_header(["HcEtac2PiPiAlg/HcEtac2PiPi.h"])
sel_m02 = Selection.new
sel_m02.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=3"   # primary pi+ and two pi+ from eta_c
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
          .kinematic_fit([:pip, :pim, :gamma, :pip, :pim]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_m02.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 35 for final states with only charged or K_S0 particles (about 80% efficiency). Loose default (200) used in BOSS; the published cut is applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-50 MeV/c^2 around the nominal eta_c mass (charged-only final state).")
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
  config.events        = 50_000
  config.decay_card    = decay_card_m03
  config.cross_section = :default
end

alg_m03 = Algorithm.new("HcEtac2KK")
alg_m03.set_header(["HcEtac2KKAlg/HcEtac2KK.h"])
sel_m03 = Selection.new
sel_m03.select_track {
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
             nkp ">=1"; nkm ">=1"
           }
          .kinematic_fit([:pip, :pim, :gamma, :kp, :km, :kp, :km]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_m03.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 35 for charged-only final states; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-50 MeV/c^2 around the nominal eta_c mass (charged-only final state).")
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
  config.sample_name   = "sig_pipihc_etac_KKpipi"
  config.events        = 50_000
  config.decay_card    = decay_card_m04
  config.cross_section = :default
end

alg_m04 = Algorithm.new("HcEtacKKPiPi")
alg_m04.set_header(["HcEtacKKPiPiAlg/HcEtacKKPiPi.h"])
sel_m04 = Selection.new
sel_m04.select_track {
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
          .kinematic_fit([:pip, :pim, :gamma, :kp, :km, :pip, :pim]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_m04.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 35 for charged-only final states; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-50 MeV/c^2 around the nominal eta_c mass (charged-only final state).")
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
  config.events        = 50_000
  config.decay_card    = decay_card_m05
  config.cross_section = :default
end

alg_m05 = Algorithm.new("HcEtacPPbarPiPi")
alg_m05.set_header(["HcEtacPPbarPiPiAlg/HcEtacPPbarPiPi.h"])
sel_m05 = Selection.new
sel_m05.select_track {
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
alg_m05.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 35 for charged-only final states; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-50 MeV/c^2 around the nominal eta_c mass (charged-only final state).")
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
  config.events        = 50_000
  config.decay_card    = decay_card_m06
  config.cross_section = :default
end

alg_m06 = Algorithm.new("HcEtac3PiPi")
alg_m06.set_header(["HcEtac3PiPiAlg/HcEtac3PiPi.h"])
sel_m06 = Selection.new
sel_m06.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=4"
             nChrn     ">=4"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=1"
           }
          .kinematic_fit([:pip, :pim, :gamma, :pip, :pim, :pip, :pim]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_m06.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 35 for charged-only final states; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-50 MeV/c^2 around the nominal eta_c mass (charged-only final state).")
alg_m06.with_decay_card(decay_card_m06).apply(sel_m06)
alg_m06.execute_on(scan_points + incMC_points + exMC_m06)

# =============================================================================
# Mode 07 : eta_c -> K+ K- 2 (pi+ pi-)
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
  config.events        = 50_000
  config.decay_card    = decay_card_m07
  config.cross_section = :default
end

alg_m07 = Algorithm.new("HcEtacKK2PiPi")
alg_m07.set_header(["HcEtacKK2PiPiAlg/HcEtacKK2PiPi.h"])
sel_m07 = Selection.new
sel_m07.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=4"
             nChrn     ">=4"
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
alg_m07.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 35 for charged-only final states; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-50 MeV/c^2 around the nominal eta_c mass (charged-only final state).")
alg_m07.with_decay_card(decay_card_m07).apply(sel_m07)
alg_m07.execute_on(scan_points + incMC_points + exMC_m07)

# =============================================================================
# Mode 08 : eta_c -> K_S0 K+ pi-
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
  config.events        = 50_000
  config.decay_card    = decay_card_m08
  config.cross_section = :default
end

alg_m08 = Algorithm.new("HcEtacKSKPi")
alg_m08.set_header(["HcEtacKSKPiAlg/HcEtacKSKPi.h"])
sel_m08 = Selection.new
sel_m08.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=3"   # primary pi+, K+, pi+ from K_S0
             nChrn     ">=3"   # primary pi-, pi- from K_S0, pi-
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
          .secondary_vertex_fit([:pip, :pim]) {
             build_virtual_particle(:K_S0).by_minimizing_mass_difference
             remove_used_particle_from_candidate_list
           }
          .kinematic_fit([:pip, :pim, :gamma, :K_S0, :kp, :pim]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_m08.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 35 for charged-only / K_S0 final states; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-50 MeV/c^2 around the nominal eta_c mass (K_S0 final state).")
      .note(:ks0_selection, "K_S0 -> pi+ pi- candidates are selected as described in Ref. [13].")
alg_m08.with_decay_card(decay_card_m08).apply(sel_m08)
alg_m08.execute_on(scan_points + incMC_points + exMC_m08)

# =============================================================================
# Mode 09 : eta_c -> K_S0 K+ pi- pi+ pi-
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
  config.sample_name   = "sig_pipihc_etac_KSKpipi"
  config.events        = 50_000
  config.decay_card    = decay_card_m09
  config.cross_section = :default
end

alg_m09 = Algorithm.new("HcEtacKSKPiPi")
alg_m09.set_header(["HcEtacKSKPiPiAlg/HcEtacKSKPiPi.h"])
sel_m09 = Selection.new
sel_m09.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=4"
             nChrn     ">=4"
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
          .secondary_vertex_fit([:pip, :pim]) {
             build_virtual_particle(:K_S0).by_minimizing_mass_difference
             remove_used_particle_from_candidate_list
           }
          .kinematic_fit([:pip, :pim, :gamma, :K_S0, :kp, :pim, :pip, :pim]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_m09.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 35 for charged-only / K_S0 final states; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-50 MeV/c^2 around the nominal eta_c mass (K_S0 final state).")
      .note(:ks0_selection, "K_S0 -> pi+ pi- candidates are selected as described in Ref. [13].")
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
  config.events        = 50_000
  config.decay_card    = decay_card_m10
  config.cross_section = :default
end

alg_m10 = Algorithm.new("HcEtacKKPi0")
alg_m10.set_header(["HcEtacKKPi0Alg/HcEtacKKPi0.h"])
sel_m10 = Selection.new
sel_m10.select_track {
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
             nGam              ">=3"   # gamma from h_c plus the two photons of pi0
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             nkp ">=1"; nkm ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 25            # 1C mass constraint; |M(gamma gamma) - m(pi0)| < 15 MeV/c^2
             npi0 ">=1"
           }
          .kinematic_fit([:pip, :pim, :gamma, :kp, :km, :pi0]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_m10.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 20 for final states containing pi0 or eta (about 70% efficiency); loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C (the 1C chi2 is the mass constraint of the two photons of each pi0/eta); resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:pi0_mass_window, "pi0 (eta) candidate: |M(gamma gamma) - m(pi0)| < 15 MeV/c^2 (|M(gamma gamma) - m(eta)| < 15 MeV/c^2).")
      .note(:eta_c_mass_window, "eta_c candidate window +/-45 MeV/c^2 around the nominal eta_c mass for final states containing pi0 or eta.")
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
  config.events        = 50_000
  config.decay_card    = decay_card_m11
  config.cross_section = :default
end

alg_m11 = Algorithm.new("HcEtacPPbarPi0")
alg_m11.set_header(["HcEtacPPbarPi0Alg/HcEtacPPbarPi0.h"])
sel_m11 = Selection.new
sel_m11.select_track {
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
             nGam              ">=3"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :proton, against: [:kaon, :pion]
             identify :pion,   against: [:kaon, :proton]
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
alg_m11.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 20 for final states containing pi0 or eta; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:pi0_mass_window, "pi0 candidate: |M(gamma gamma) - m(pi0)| < 15 MeV/c^2.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-45 MeV/c^2 around the nominal eta_c mass.")
alg_m11.with_decay_card(decay_card_m11).apply(sel_m11)
alg_m11.execute_on(scan_points + incMC_points + exMC_m11)

# =============================================================================
# Mode 12 : eta_c -> pi+ pi- eta
# =============================================================================
decay_card_m12 = <<~DECAYCARD
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

exMC_m12 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_pipieta"
  config.events        = 50_000
  config.decay_card    = decay_card_m12
  config.cross_section = :default
end

alg_m12 = Algorithm.new("HcEtacPiPiEta")
alg_m12.set_header(["HcEtacPiPiEtaAlg/HcEtacPiPiEta.h"])
sel_m12 = Selection.new
sel_m12.select_track {
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
             nGam              ">=3"   # gamma from h_c plus the two photons of eta
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25            # |M(gamma gamma) - m(eta)| < 15 MeV/c^2
             neta ">=1"
           }
          .kinematic_fit([:pip, :pim, :gamma, :pip, :pim, :eta]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_m12.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 20 for final states containing pi0 or eta; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_mass_window, "eta candidate: |M(gamma gamma) - m(eta)| < 15 MeV/c^2.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-45 MeV/c^2 around the nominal eta_c mass.")
alg_m12.with_decay_card(decay_card_m12).apply(sel_m12)
alg_m12.execute_on(scan_points + incMC_points + exMC_m12)

# =============================================================================
# Mode 13 : eta_c -> K+ K- eta
# =============================================================================
decay_card_m13 = <<~DECAYCARD
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

exMC_m13 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_KKeta"
  config.events        = 50_000
  config.decay_card    = decay_card_m13
  config.cross_section = :default
end

alg_m13 = Algorithm.new("HcEtacKKEta")
alg_m13.set_header(["HcEtacKKEtaAlg/HcEtacKKEta.h"])
sel_m13 = Selection.new
sel_m13.select_track {
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
             nGam              ">=3"
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
alg_m13.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 20 for final states containing pi0 or eta; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_mass_window, "eta candidate: |M(gamma gamma) - m(eta)| < 15 MeV/c^2.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-45 MeV/c^2 around the nominal eta_c mass.")
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
  config.events        = 50_000
  config.decay_card    = decay_card_m14
  config.cross_section = :default
end

alg_m14 = Algorithm.new("HcEtac2PiPiEta")
alg_m14.set_header(["HcEtac2PiPiEtaAlg/HcEtac2PiPiEta.h"])
sel_m14 = Selection.new
sel_m14.select_track {
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
             nGam              ">=3"
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
alg_m14.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 20 for final states containing pi0 or eta; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:eta_mass_window, "eta candidate: |M(gamma gamma) - m(eta)| < 15 MeV/c^2.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-45 MeV/c^2 around the nominal eta_c mass.")
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
  config.events        = 50_000
  config.decay_card    = decay_card_m15
  config.cross_section = :default
end

alg_m15 = Algorithm.new("HcEtacPiPi2Pi0")
alg_m15.set_header(["HcEtacPiPi2Pi0Alg/HcEtacPiPi2Pi0.h"])
sel_m15 = Selection.new
sel_m15.select_track {
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
             nGam              ">=5"   # gamma from h_c plus four photons of the two pi0
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 25
             npi0 ">=2"
           }
          .kinematic_fit([:pip, :pim, :gamma, :pi0, :pi0]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_m15.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 20 for final states containing pi0 or eta; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:pi0_mass_window, "pi0 candidates: |M(gamma gamma) - m(pi0)| < 15 MeV/c^2.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-45 MeV/c^2 around the nominal eta_c mass.")
alg_m15.with_decay_card(decay_card_m15).apply(sel_m15)
alg_m15.execute_on(scan_points + incMC_points + exMC_m15)

# =============================================================================
# Mode 16 : eta_c -> 2 (pi+ pi-) pi0 pi0
# =============================================================================
decay_card_m16 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi+ pi- pi0 pi0   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m16 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_pipihc_etac_2pipi2pi0"
  config.events        = 50_000
  config.decay_card    = decay_card_m16
  config.cross_section = :default
end

alg_m16 = Algorithm.new("HcEtac2PiPi2Pi0")
alg_m16.set_header(["HcEtac2PiPi2Pi0Alg/HcEtac2PiPi2Pi0.h"])
sel_m16 = Selection.new
sel_m16.select_track {
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
             nGam              ">=5"
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
alg_m16.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 20 for final states containing pi0 or eta; loose default (200) used in BOSS, published cut applied in ROOT.")
      .note(:combination_selection, "Combination minimising chi2 = chi2_4C + sum_i chi2_PID(i) + chi2_1C; resolved in BOSS by the minimum chi2_4C.")
      .note(:recoil_mass_requirements, "M_recoil(pi+pi-) in [3.45, 3.65] GeV/c^2 and M_recoil(gamma pi+pi-) in [2.8, 3.2] GeV/c^2 are required for at least one combination.")
      .note(:pi0_mass_window, "pi0 candidates: |M(gamma gamma) - m(pi0)| < 15 MeV/c^2.")
      .note(:eta_c_mass_window, "eta_c candidate window +/-45 MeV/c^2 around the nominal eta_c mass.")
alg_m16.with_decay_card(decay_card_m16).apply(sel_m16)
alg_m16.execute_on(scan_points + incMC_points + exMC_m16)
