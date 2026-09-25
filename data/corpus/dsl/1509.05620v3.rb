### Dataset description ###
data_4226 = DatasetManager.real_data.find("703_4230")   # sqrt(s)=4.226 GeV (1092/pb)
data_4257 = DatasetManager.real_data.find("703_4260")   # sqrt(s)=4.257 GeV (826/pb)
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4257 = DatasetManager.inclusive_mc.find("703_4260")

# ---------------- Decay cards ----------------
# Signal Mode I: e+e- -> D+ D*- pi0, D*- -> Dbar0 pi- ; D+ / Dbar0 in various sub-modes
decay_card_ModeI = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D+ D*- pi0                     PHSP;
  Enddecay

  Decay D*-
  1.0000 anti-D0 pi-                    VSS;
  Enddecay

  Decay D+
  0.2000 K- pi+ pi+                     PHSP;
  0.2000 K- pi+ pi+ pi0                 PHSP;
  0.2000 K_S0 pi+                       PHSP;
  0.2000 K_S0 pi+ pi0                   PHSP;
  0.2000 K_S0 pi+ pi+ pi-               PHSP;
  Enddecay

  Decay anti-D0
  0.3333 K+ pi-                         PHSP;
  0.3333 K+ pi- pi0                     PHSP;
  0.3334 K+ pi- pi- pi+                 PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                        PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                    PHSP;
  Enddecay

  End
DECAYCARD

# Signal Mode II: e+e- -> D0 Dbar*0 pi0, Dbar*0 -> Dbar0 pi0
decay_card_ModeII = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D0 anti-D*0 pi0                PHSP;
  Enddecay

  Decay anti-D*0
  1.0000 anti-D0 pi0                    VSS;
  Enddecay

  Decay D0
  0.3333 K- pi+                         PHSP;
  0.3333 K- pi+ pi0                     PHSP;
  0.3334 K- pi+ pi+ pi-                 PHSP;
  Enddecay

  Decay anti-D0
  0.3333 K+ pi-                         PHSP;
  0.3333 K+ pi- pi0                     PHSP;
  0.3334 K+ pi- pi- pi+                 PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                        PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                    PHSP;
  Enddecay

  End
DECAYCARD

# ---------------- Exclusive MC ----------------
exMC_ModeI = DatasetManager.create_exclusive_mc_for([data_4226, data_4257]) do |config|
  config.sample_name = "Zc3885_neutral_DpDstm_pi0"
  config.events = 200000
  config.decay_card = decay_card_ModeI
  config.cross_section = :default
end
exMC_ModeI.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

exMC_ModeII = DatasetManager.create_exclusive_mc_for([data_4226, data_4257]) do |config|
  config.sample_name = "Zc3885_neutral_D0Dst0bar_pi0"
  config.events = 200000
  config.decay_card = decay_card_ModeII
  config.cross_section = :default
end
exMC_ModeII.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

# ================= Event Selection: Mode I (D+ D*- pi0) =================
alg_ModeI = Algorithm.new("Zc3885Neutral_ModeI")
alg_ModeI.set_header(["Zc3885Neutral_ModeIAlg/Zc3885Neutral_ModeI.h"])
         .set_constant({"ECMS" => [:double, 4.226]})

sel_ModeI = Selection.new
sel_ModeI.select_track {
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     ">=2"
            nChrn     ">=1"
          }
         .select_photon {
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            tdc_emc_start     0
            tdc_emc_end       14
            nGam              ">=2"
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion, :proton]
            identify :pion, against: [:kaon, :proton]
          }
          # Reconstruct primary pi0 (photons not used to form D mesons; |M(gg)-m(pi0)| window)
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 200
            npi0     ">=1"
          }
          # Reconstruct D+ via one of five sub-modes; mass window (1.840,1.880) GeV, then mass-constrained KF
         .kinematic_fit([:km, :pip, :pip, :kp, :pim, :pi0]) {
            invariant_mass_of(:km, :pip, :pip).within(1.840, 1.880)
            invariant_mass_of(:kp, :pim, :pi0).within(1.840, 1.880)
            invariant_mass_of(:km, :pip, :pip).constrain_to_nominal_mass_of(:"D+")
            invariant_mass_of(:kp, :pim, :pi0).constrain_to_nominal_mass_of(:"anti-D0")
            chi2_cut 100
          }
          # Primary pi0 2C-fit and (D pi0) kinematic requirements: nominal KF constrains
          # M(gg)->m(pi0) and RM(pi0 D Dbar)->m(pi); RM(D pi0) taken as :D*-
         .kinematic_fit([:"D+", :"anti-D0", :gamma, :gamma]) {
            nominal
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            invariant_mass_of(:"D+", :"anti-D0", :gamma, :gamma).recoil_mass.constrain_to_nominal_mass_of(:"pi-")
            chi2_cut 200
          }

alg_ModeI.note(:soft_pion_missing,
               "The soft pion from D*- -> Dbar0 pi- is not required in reconstruction; " \
               "detection is based on the D+ Dbar0 pair.")
         .note(:d_common_vertex_fit,
               "Charged tracks in each D candidate (excluding those from K_S0 -> pi+pi-) " \
               "are constrained to a common vertex with chi2 < 100.")
         .note(:pi0_window,
               "M(gamma gamma) required in (0.120, 0.150) GeV before the 2C KF.")
         .note(:multi_candidate_selection,
               "Among multiple DD combinations, keep the candidate with minimum sum " \
               "of chi2_D + chi2_Dbar from the mass-constrained fits.")
         .note(:Dpi0_mass_cut,
               "Require M(D+ pi0) > 2.1 GeV and M(Dbar0 pi0) > 2.1 GeV.")
         .note(:RM_Dpi0_window,
               "Require |RM(D pi0) - m(D*)| < 36 MeV/c^2 to select (D Dbar*)^0 pi0 candidates.")
         .with_decay_card(decay_card_ModeI)
         .apply(sel_ModeI)

alg_ModeI.execute_on([data_4226, data_4257, incMC_4226, incMC_4257] + exMC_ModeI)

# ================= Event Selection: Mode II (D0 Dbar*0 pi0, Dbar*0->Dbar0 pi0) =================
alg_ModeII = Algorithm.new("Zc3885Neutral_ModeII")
alg_ModeII.set_header(["Zc3885Neutral_ModeIIAlg/Zc3885Neutral_ModeII.h"])
          .set_constant({"ECMS" => [:double, 4.226]})

sel_ModeII = Selection.new
sel_ModeII.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=1"
             nChrn     ">=1"
           }
          .select_photon {
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             tdc_emc_start     0
             tdc_emc_end       14
             nGam              ">=2"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             identify :pion, against: [:kaon, :proton]
           }
           # Primary pi0 reconstruction (photons not used to form D mesons)
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 200
             npi0     ">=1"
           }
           # Reconstruct D0 and Dbar0, mass window (1.840, 1.880) GeV; mass-constrained KF
          .kinematic_fit([:km, :pip, :kp, :pim]) {
             invariant_mass_of(:km, :pip).within(1.840, 1.880)
             invariant_mass_of(:kp, :pim).within(1.840, 1.880)
             invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)
             invariant_mass_of(:kp, :pim).constrain_to_nominal_mass_of(:"anti-D0")
             chi2_cut 100
           }
           # Primary pi0 2C fit: M(gg)=m(pi0) and RM(pi0 D D)=m(pi)
          .kinematic_fit([:D0, :"anti-D0", :gamma, :gamma]) {
             nominal
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             invariant_mass_of(:D0, :"anti-D0", :gamma, :gamma).recoil_mass.constrain_to_nominal_mass_of(:pi0)
             chi2_cut 60
           }
           # Competing hypothesis: Dbar*0 -> Dbar0 gamma; constrain RM(pi0 D0 Dbar0)=0
          .kinematic_fit([:D0, :"anti-D0", :gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             invariant_mass_of(:D0, :"anti-D0", :gamma, :gamma).recoil_mass.within(-0.001, 0.001)
           }

alg_ModeII.note(:pi0_window,
                "M(gamma gamma) required in (0.120, 0.150) GeV before the 2C KF.")
          .note(:d_common_vertex_fit,
                "Charged tracks in each D candidate constrained to common vertex with chi2 < 100.")
          .note(:multi_candidate_selection,
                "Best DD combination taken as minimum sum of chi2_D + chi2_Dbar.")
          .note(:radiative_veto,
                "Suppress Dbar*0 -> Dbar0 gamma background by requiring chi2_2C(pi0) < 60 " \
                "AND chi2_2C(gamma) > 20, where chi2_2C(gamma) is the 2C KF chi2 for the " \
                "photon hypothesis (RM(pi0 D0 Dbar0) = 0).")
          .note(:Dpi0_mass_cut,
                "Require M(D0 pi0) > 2.1 GeV for both D0 and Dbar0.")
          .note(:RM_Dpi0_window,
                "Require |RM(D pi0) - m(D*)| < 36 MeV/c^2.")
          .with_decay_card(decay_card_ModeII)
          .apply(sel_ModeII)

alg_ModeII.execute_on([data_4226, data_4257, incMC_4226, incMC_4257] + exMC_ModeII)
