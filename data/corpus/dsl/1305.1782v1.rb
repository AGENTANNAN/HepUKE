# =============================================================================
# BESIII : Search for psi(3770) and psi(4040) decays to baryonic final states
# arXiv:1305.1782v1
#   Lambda Lambdabar pi+pi- , Lambda Lambdabar pi0 , Lambda Lambdabar eta ,
#   Sigma+ Sigmabar- , Sigma0 Sigmabar0 , Xi- Xibar+ , Xi0 Xibar0
# Data : 2.9 fb^-1 @ 3.773 GeV, 482 pb^-1 @ 4.009 GeV,
#        23 pb^-1 @ 3.542/3.554/3.561/3.600 GeV, 44 pb^-1 @ 3.650 GeV
# BOSS part only : dataset preparation + event selection up to the 4C fit.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets
### ---------------------------------------------------------------------------
psi3770_data  = DatasetManager.real_data.find("712_3773")   # psi(3770), 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")
psi4040_data  = DatasetManager.real_data.find("703_4009")   # psi(4040), 4.009 GeV
psi4040_incMC = DatasetManager.inclusive_mc.find("703_4009")
# Non-resonant continuum data used for the qqbar background subtraction
cont_data     = DatasetManager.real_data.find("712_3554")   # 3.554 GeV (23 pb-1 @ 3.542-3.600)
cont3650_data = DatasetManager.real_data.find("709_3650")   # 3.650 GeV (44 pb-1)

common_datasets = [psi3770_data, psi4040_data, psi3770_incMC, psi4040_incMC,
                   cont_data, cont3650_data]

# =============================================================================
# Mode 1 : psi(3770) -> Lambda Lambdabar pi+ pi-
# =============================================================================
decay_card_mode1 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 pi+ pi-   PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-   HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+   HypWK;
  Enddecay

  End
DECAYCARD

exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psi3770_LLpipi"
  config.related_dataset = psi3770_data
  config.events          = 50_000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

alg_mode1 = Algorithm.new("PsipBaryonicLLpipi")
alg_mode1.set_header(["PsipBaryonicLLpipiAlg/PsipBaryonicLLpipi.h"])
         .set_constant({"ECMS" => [:double, 3.773]})

sel_mode1 = Selection.new
sel_mode1.select_track {
            cos_theta 0.93      # |cos(theta)| < 0.93 w.r.t. the e+ direction
            Vz        100.0     # no primary-vertex requirement on the Lambda daughters
            Vr        10.0
            nChrp     ">=3"     # 3 positive tracks : p, pi+ (from Lambdabar), pi+
            nChrn     ">=3"     # 3 negative tracks : pbar, pi- (from Lambda), pi-
            nNet      "==0"
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]   # CL_p > 0.001, > CL_pi, > CL_K
            identify :pion,   against: [:kaon, :proton] # CL_pi > CL_K
            nprp ">=1"; nprm ">=1"
            npip ">=1"; npim ">=1"
          }
         # p / pbar transverse momentum must exceed 300 MeV/c
         .for_each(:prp) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .for_each(:prm) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         # Lambda -> p pi- secondary vertex fit
         .secondary_vertex_fit([:prp, :pim]) {
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         # Lambdabar -> pbar pi+ secondary vertex fit
         .secondary_vertex_fit([:prm, :pip]) {
            build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         # 4C kinematic fit over the full final state
         .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_mode1.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 60 (loose default 200 used in BOSS; the published tight cut is applied in the ROOT stage).")
         .note(:secondary_vertex_fit_algorithm, "For Lambda (and Xi- -> Lambda pi-) the secondary vertex fit additionally imposes the kinematic constraint between production and decay vertex using the run-by-run averaged interaction point; only the p pi- vertex fit is used for the baryon-pair modes Sigma0 Sigmabar0, Xi- Xibar+, Xi0 Xibar0.")
         .note(:intermediate_mass_windows, "Loose requirement |M(p pi-) - M(Lambda)| < 40 MeV/c^2 for all modes containing a Lambda. Signal windows (about 3 sigma): Lambda 1.107-1.124, pi0 115-150 MeV, eta 515-569 MeV, Sigma+ 1.164-1.206, Sigma0 1.178-1.205, Xi- 1.305-1.337, Xi0 1.281-1.330 GeV/c^2. Sidebands (5-8 sigma) and the 2D signal boxes are applied in ROOT.")
         .note(:no_primary_vertex_requirement, "Tracks used to reconstruct Lambda, Sigma+, Sigma0, Xi- and Xi0 are not required to satisfy a primary-vertex requirement.")
         .note(:signal_mc_energy_points, "Signal MC is also generated at sqrt(s)=4.009 GeV (psi(4040) as parent) and for the continuum points 3.542/3.554/3.561/3.600/3.650 GeV, 50000 events per channel.")

alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)
alg_mode1.execute_on(common_datasets + [exMC_mode1])

# =============================================================================
# Mode 2 : psi(3770) -> Lambda Lambdabar pi0
# =============================================================================
decay_card_mode2 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 pi0   PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-   HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+   HypWK;
  Enddecay

  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay

  End
DECAYCARD

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psi3770_LLpi0"
  config.related_dataset = psi3770_data
  config.events          = 50_000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

alg_mode2 = Algorithm.new("PsipBaryonicLLpi0")
alg_mode2.set_header(["PsipBaryonicLLpi0Alg/PsipBaryonicLLpi0.h"])
         .set_constant({"ECMS" => [:double, 3.773]})

sel_mode2 = Selection.new
sel_mode2.select_track {
            cos_theta 0.93
            Vz        100.0
            Vr        10.0
            nChrp     ">=1"
            nChrn     ">=1"
            nNet      "==0"
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025   # barrel  (|cos theta| < 0.8)
            energyThreshold_e 0.050   # end caps (0.86 < |cos theta| < 0.92)
            angle_to_track    10.0
            nGam              ">=2"
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]
            identify :pion,   against: [:kaon, :proton]
            nprp ">=1"; nprm ">=1"
            npip ">=1"; npim ">=1"
          }
         .for_each(:prp) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .for_each(:prm) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .secondary_vertex_fit([:prp, :pim]) {
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         .secondary_vertex_fit([:prm, :pip]) {
            build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         # pi0 -> gamma gamma
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
          }
         .kinematic_fit([:Lambda, :Lambda_bar, :pi0]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_mode2.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 60 (loose default 200 used in BOSS; the published tight cut is applied in the ROOT stage).")
         .note(:secondary_vertex_fit_algorithm, "Lambda reconstructed with the vertex fit of p and pi- plus the secondary vertex fit algorithm constraining the production and decay vertices with the run-by-run averaged interaction point.")
         .note(:intermediate_mass_windows, "|M(p pi-) - M(Lambda)| < 40 MeV/c^2; pi0 signal window 115-150 MeV, Lambda signal window 1.107-1.124 GeV/c^2. Normalised pi0 sideband events are subtracted from the signal region.")
         .note(:no_primary_vertex_requirement, "Tracks used to reconstruct Lambda are not required to satisfy a primary-vertex requirement.")
         .note(:signal_mc_energy_points, "Signal MC is also generated at sqrt(s)=4.009 GeV (psi(4040) as parent) and for the continuum points 3.542/3.554/3.561/3.600/3.650 GeV, 50000 events per channel.")

alg_mode2.with_decay_card(decay_card_mode2).apply(sel_mode2)
alg_mode2.execute_on(common_datasets + [exMC_mode2])

# =============================================================================
# Mode 3 : psi(3770) -> Lambda Lambdabar eta
# =============================================================================
decay_card_mode3 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 eta   PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-   HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+   HypWK;
  Enddecay

  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay

  End
DECAYCARD

exMC_mode3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psi3770_LLeta"
  config.related_dataset = psi3770_data
  config.events          = 50_000
  config.decay_card      = decay_card_mode3
  config.cross_section   = :default
end

alg_mode3 = Algorithm.new("PsipBaryonicLLeta")
alg_mode3.set_header(["PsipBaryonicLLetaAlg/PsipBaryonicLLeta.h"])
         .set_constant({"ECMS" => [:double, 3.773]})

sel_mode3 = Selection.new
sel_mode3.select_track {
            cos_theta 0.93
            Vz        100.0
            Vr        10.0
            nChrp     ">=1"
            nChrn     ">=1"
            nNet      "==0"
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track    10.0
            nGam              ">=2"
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]
            identify :pion,   against: [:kaon, :proton]
            nprp ">=1"; nprm ">=1"
            npip ">=1"; npim ">=1"
          }
         .for_each(:prp) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .for_each(:prm) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .secondary_vertex_fit([:prp, :pim]) {
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         .secondary_vertex_fit([:prm, :pip]) {
            build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         # eta -> gamma gamma
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
            chi2_cut 25
            neta ">=1"
          }
         .kinematic_fit([:Lambda, :Lambda_bar, :eta]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_mode3.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 60 (loose default 200 used in BOSS; the published tight cut is applied in the ROOT stage).")
         .note(:secondary_vertex_fit_algorithm, "Lambda reconstructed with the vertex fit of p and pi- plus the secondary vertex fit algorithm constraining the production and decay vertices with the run-by-run averaged interaction point.")
         .note(:intermediate_mass_windows, "|M(p pi-) - M(Lambda)| < 40 MeV/c^2; eta signal window 515-569 MeV, Lambda signal window 1.107-1.124 GeV/c^2. Normalised eta sideband events are subtracted from the signal region.")
         .note(:no_primary_vertex_requirement, "Tracks used to reconstruct Lambda are not required to satisfy a primary-vertex requirement.")
         .note(:signal_mc_energy_points, "Signal MC is also generated at sqrt(s)=4.009 GeV (psi(4040) as parent) and for the continuum points 3.542/3.554/3.561/3.600/3.650 GeV, 50000 events per channel.")

alg_mode3.with_decay_card(decay_card_mode3).apply(sel_mode3)
alg_mode3.execute_on(common_datasets + [exMC_mode3])

# =============================================================================
# Mode 4 : psi(3770) -> Sigma+ Sigmabar-
# =============================================================================
decay_card_mode4 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Sigma+ anti-Sigma-   PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0   PHSP;
  Enddecay

  Decay anti-Sigma-
  1.0000 anti-p- pi0   PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay

  End
DECAYCARD

exMC_mode4 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psi3770_SS"
  config.related_dataset = psi3770_data
  config.events          = 50_000
  config.decay_card      = decay_card_mode4
  config.cross_section   = :default
end

alg_mode4 = Algorithm.new("PsipBaryonicSS")
alg_mode4.set_header(["PsipBaryonicSSAlg/PsipBaryonicSS.h"])
         .set_constant({"ECMS" => [:double, 3.773]})

sel_mode4 = Selection.new
sel_mode4.select_track {
            cos_theta 0.93
            Vz        100.0
            Vr        10.0
            nChrp     ">=1"
            nChrn     ">=1"
            nNet      "==0"
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track    10.0
            nGam              ">=4"   # Sigma+ -> p pi0 and Sigmabar- -> pbar pi0
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]
            identify :pion,   against: [:kaon, :proton]
            nprp ">=1"; nprm ">=1"
          }
         .for_each(:prp) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .for_each(:prm) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         # two independent pi0 -> gamma gamma
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=2"
          }
         .kinematic_fit([:prp, :prm, :pi0, :pi0]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_mode4.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 60 (loose default 200 used in BOSS; the published tight cut is applied in the ROOT stage).")
         .note(:multiple_combination_selection, "For the two-pi0 final states the candidate minimising R(pi0) = sqrt[(M(gamma gamma)_1 - M(pi0))^2 + (M(gamma gamma)_2 - M(pi0))^2] is retained; for the baryon-pair modes formed from p pbar pi0 pi0 / Lambda Lambdabar gamma gamma / Lambda Lambdabar pi0 pi0 the minimum of R(j) = sqrt[(M(i)-M(j))^2 + (M(i')-M(j))^2] is retained. In the DSL the combination is resolved by the minimum chi2_4C of the kinematic fit; the R-based choice is made in ROOT.")
         .note(:intermediate_mass_windows, "pi0 window 115-150 MeV, Sigma+ window 1.164-1.206 GeV/c^2 (M(p gamma gamma)), Sigma0 window 1.178-1.205 GeV/c^2.")
         .note(:no_primary_vertex_requirement, "Tracks used to reconstruct Sigma+ / Sigma0 are not required to satisfy a primary-vertex requirement.")
         .note(:signal_mc_energy_points, "Signal MC is also generated at sqrt(s)=4.009 GeV (psi(4040) as parent) and for the continuum points 3.542/3.554/3.561/3.600/3.650 GeV, 50000 events per channel.")

alg_mode4.with_decay_card(decay_card_mode4).apply(sel_mode4)
alg_mode4.execute_on(common_datasets + [exMC_mode4])

# =============================================================================
# Mode 5 : psi(3770) -> Sigma0 Sigmabar0  (Sigma0 -> Lambda gamma)
# =============================================================================
decay_card_mode5 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Sigma0 anti-Sigma0   PHSP;
  Enddecay

  Decay Sigma0
  1.0000 Lambda0 gamma   PHSP;
  Enddecay

  Decay anti-Sigma0
  1.0000 anti-Lambda0 gamma   PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-   HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+   HypWK;
  Enddecay

  End
DECAYCARD

exMC_mode5 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psi3770_S0S0"
  config.related_dataset = psi3770_data
  config.events          = 50_000
  config.decay_card      = decay_card_mode5
  config.cross_section   = :default
end

alg_mode5 = Algorithm.new("PsipBaryonicS0S0")
alg_mode5.set_header(["PsipBaryonicS0S0Alg/PsipBaryonicS0S0.h"])
         .set_constant({"ECMS" => [:double, 3.773]})

sel_mode5 = Selection.new
sel_mode5.select_track {
            cos_theta 0.93
            Vz        100.0
            Vr        10.0
            nChrp     ">=1"
            nChrn     ">=1"
            nNet      "==0"
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track    10.0
            nGam              ">=2"
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]
            identify :pion,   against: [:kaon, :proton]
            nprp ">=1"; nprm ">=1"
            npip ">=1"; npim ">=1"
          }
         .for_each(:prp) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .for_each(:prm) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         # only the p pi- vertex fit is used for the baryon-pair modes
         .secondary_vertex_fit([:prp, :pim]) {
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         .secondary_vertex_fit([:prm, :pip]) {
            build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_mode5.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 60 (loose default 200 used in BOSS; the published tight cut is applied in the ROOT stage).")
         .note(:secondary_vertex_fit_algorithm, "For this baryon-pair mode only the vertex fit of p and pi- is used to reconstruct the Lambda (no production-vertex constraint).")
         .note(:multiple_combination_selection, "The two Lambdas are formed from p pbar pi+ pi- and two photons; the minimum of R(Sigma0) = sqrt[(M(Lambda gamma)-M(Sigma0))^2 + (M(Lambdabar gamma)-M(Sigma0))^2] selects the assignment. Resolved here by the minimum chi2_4C of the kinematic fit.")
         .note(:intermediate_mass_windows, "|M(p pi-) - M(Lambda)| < 40 MeV/c^2; Lambda 1.107-1.124 GeV/c^2, Sigma0 1.178-1.205 GeV/c^2 (M(p pi- gamma)).")
         .note(:no_primary_vertex_requirement, "Tracks used to reconstruct Lambda and Sigma0 are not required to satisfy a primary-vertex requirement.")
         .note(:signal_mc_energy_points, "Signal MC is also generated at sqrt(s)=4.009 GeV (psi(4040) as parent) and for the continuum points 3.542/3.554/3.561/3.600/3.650 GeV, 50000 events per channel.")

alg_mode5.with_decay_card(decay_card_mode5).apply(sel_mode5)
alg_mode5.execute_on(common_datasets + [exMC_mode5])

# =============================================================================
# Mode 6 : psi(3770) -> Xi- Xibar+  (Xi- -> Lambda pi-)
# =============================================================================
decay_card_mode6 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Xi- anti-Xi+   PHSP;
  Enddecay

  Decay Xi-
  1.0000 Lambda0 pi-   PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+   PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-   HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+   HypWK;
  Enddecay

  End
DECAYCARD

exMC_mode6 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psi3770_XiXi"
  config.related_dataset = psi3770_data
  config.events          = 50_000
  config.decay_card      = decay_card_mode6
  config.cross_section   = :default
end

alg_mode6 = Algorithm.new("PsipBaryonicXiXi")
alg_mode6.set_header(["PsipBaryonicXiXiAlg/PsipBaryonicXiXi.h"])
         .set_constant({"ECMS" => [:double, 3.773]})

sel_mode6 = Selection.new
sel_mode6.select_track {
            cos_theta 0.93
            Vz        100.0
            Vr        10.0
            nChrp     ">=3"   # p, pi+ (Lambdabar), pi+ (anti-Xi+)
            nChrn     ">=3"   # pbar, pi- (Lambda), pi- (Xi-)
            nNet      "==0"
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]
            identify :pion,   against: [:kaon, :proton]
            nprp ">=1"; nprm ">=1"
            npip ">=1"; npim ">=1"
          }
         .for_each(:prp) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .for_each(:prm) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .secondary_vertex_fit([:prp, :pim]) {
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         .secondary_vertex_fit([:prm, :pip]) {
            build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         # remaining pi- and pi+ are the daughters of Xi- and anti-Xi+
         .kinematic_fit([:Lambda, :pim, :Lambda_bar, :pip]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_mode6.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 60 (loose default 200 used in BOSS; the published tight cut is applied in the ROOT stage).")
         .note(:secondary_vertex_fit_algorithm, "Only the p pi- vertex fit reconstructs the Lambda for this baryon-pair mode; additionally the Lambda and the pi- are fitted to a common vertex to form the Xi- (and likewise for the anti-Xi+).")
         .note(:intermediate_mass_windows, "|M(p pi-) - M(Lambda)| < 40 MeV/c^2; Lambda 1.107-1.124 GeV/c^2, Xi- 1.305-1.337 GeV/c^2 (M(p pi- pi-)).")
         .note(:no_primary_vertex_requirement, "Tracks used to reconstruct Lambda and Xi- are not required to satisfy a primary-vertex requirement.")
         .note(:signal_mc_energy_points, "Signal MC is also generated at sqrt(s)=4.009 GeV (psi(4040) as parent) and for the continuum points 3.542/3.554/3.561/3.600/3.650 GeV, 50000 events per channel.")

alg_mode6.with_decay_card(decay_card_mode6).apply(sel_mode6)
alg_mode6.execute_on(common_datasets + [exMC_mode6])

# =============================================================================
# Mode 7 : psi(3770) -> Xi0 Xibar0  (Xi0 -> Lambda pi0)
# =============================================================================
decay_card_mode7 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Xi0 anti-Xi0   PHSP;
  Enddecay

  Decay Xi0
  1.0000 Lambda0 pi0   PHSP;
  Enddecay

  Decay anti-Xi0
  1.0000 anti-Lambda0 pi0   PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-   HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+   HypWK;
  Enddecay

  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay

  End
DECAYCARD

exMC_mode7 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psi3770_Xi0Xi0"
  config.related_dataset = psi3770_data
  config.events          = 50_000
  config.decay_card      = decay_card_mode7
  config.cross_section   = :default
end

alg_mode7 = Algorithm.new("PsipBaryonicXi0Xi0")
alg_mode7.set_header(["PsipBaryonicXi0Xi0Alg/PsipBaryonicXi0Xi0.h"])
         .set_constant({"ECMS" => [:double, 3.773]})

sel_mode7 = Selection.new
sel_mode7.select_track {
            cos_theta 0.93
            Vz        100.0
            Vr        10.0
            nChrp     ">=1"
            nChrn     ">=1"
            nNet      "==0"
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track    10.0
            nGam              ">=4"   # Xi0 -> Lambda pi0 and anti-Xi0 -> Lambdabar pi0
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]
            identify :pion,   against: [:kaon, :proton]
            nprp ">=1"; nprm ">=1"
            npip ">=1"; npim ">=1"
          }
         .for_each(:prp) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .for_each(:prm) { where { (px.__pow__(2) + py.__pow__(2)) < 0.09 }; remove }
         .secondary_vertex_fit([:prp, :pim]) {
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         .secondary_vertex_fit([:prm, :pip]) {
            build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
         # two pi0 -> gamma gamma ; looser window (110-150 MeV) is applied in ROOT
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=2"
          }
         .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :pi0]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_mode7.note(:kinematic_fit_chi2, "Paper requires chi2_4C < 60 (loose default 200 used in BOSS; the published tight cut is applied in the ROOT stage).")
         .note(:secondary_vertex_fit_algorithm, "Only the p pi- vertex fit is used to reconstruct the Lambda for this baryon-pair mode.")
         .note(:multiple_combination_selection, "The baryon pair is built from Lambda Lambdabar pi0 pi0; the minimum of R(Xi0) = sqrt[(M(Lambda pi0)-M(Xi0))^2 + (M(Lambdabar pi0)-M(Xi0))^2] selects the assignment (resolved here by the minimum chi2_4C of the kinematic fit). The pi0 selection uses the looser window 110-150 MeV because of the clean signal.")
         .note(:intermediate_mass_windows, "|M(p pi-) - M(Lambda)| < 40 MeV/c^2; Lambda 1.107-1.124 GeV/c^2, Xi0 1.281-1.330 GeV/c^2 (M(p pi- gamma gamma)).")
         .note(:no_primary_vertex_requirement, "Tracks used to reconstruct Lambda and Xi0 are not required to satisfy a primary-vertex requirement.")
         .note(:signal_mc_energy_points, "Signal MC is also generated at sqrt(s)=4.009 GeV (psi(4040) as parent) and for the continuum points 3.542/3.554/3.561/3.600/3.650 GeV, 50000 events per channel.")

alg_mode7.with_decay_card(decay_card_mode7).apply(sel_mode7)
alg_mode7.execute_on(common_datasets + [exMC_mode7])
