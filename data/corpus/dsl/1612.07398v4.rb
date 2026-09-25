### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for Mode I: chi_c2 -> K+ K- pi0
decay_card_kkpi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2                    HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay chi_c2
    1.0000 K+ K- pi0                       PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                     PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for Mode II: chi_c2 -> KS K+- pi-+
decay_card_kskpi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2                    HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay chi_c2
    0.5000 K_S0 K+ pi-                     PHSP;
    0.5000 K_S0 K- pi+                     PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                         PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for Mode III: chi_c2 -> pi+ pi- pi0
decay_card_3pi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2                    HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay chi_c2
    1.0000 pi+ pi- pi0                     PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                     PHSP;
    Enddecay

    End
DECAYCARD

exMC_kkpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "chic2_KKpi0_exclusive_mc"
  config.related_dataset = psip_data
  config.events         = 200000
  config.decay_card     = decay_card_kkpi0
  config.cross_section  = :default
end

exMC_kskpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "chic2_KsKpi_exclusive_mc"
  config.related_dataset = psip_data
  config.events         = 200000
  config.decay_card     = decay_card_kskpi
  config.cross_section  = :default
end

exMC_3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "chic2_pipipi0_exclusive_mc"
  config.related_dataset = psip_data
  config.events         = 200000
  config.decay_card     = decay_card_3pi
  config.cross_section  = :default
end

########################################################################
# Mode I: psi(3686) -> gamma K+ K- pi0
########################################################################
alg_kkpi0 = Algorithm.new("Chic2ToKKpi0")
alg_kkpi0.set_header(["Chic2ToKKpi0Alg/Chic2ToKKpi0.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_kkpi0 = Selection.new
sel_kkpi0.select_track {
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     "==1"
            nChrn     "==1"
            nNet      "==0"
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=3"
         }
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
         }
         .assign({:chrgp => :kp, :chrgn => :km})
         .kinematic_fit([:gamma, :kp, :km, :gamma, :gamma]) {
            nominal
            constrain_four_momentum
            chi2_cut 80
            invariant_mass_of(:kp, :km).within(0.0, 3.0)
         }

alg_kkpi0.note(:pi0_selection,
               "pi0 candidate reconstructed from two selected photons whose invariant mass is closest to pi0 nominal mass and satisfies |M_gg - M_pi0| < 10 MeV/c^2")
         .note(:chic2_signal_window,
               "|M(K+K-pi0) - M(chi_c2)| <= 15 MeV/c^2 signal region defined offline in ROOT")
         .with_decay_card(decay_card_kkpi0)
         .apply(sel_kkpi0)

alg_kkpi0.execute_on([psip_data, psip_incMC, exMC_kkpi0])

########################################################################
# Mode II: psi(3686) -> gamma KS K+- pi-+
########################################################################
alg_kskpi = Algorithm.new("Chic2ToKsKpi")
alg_kskpi.set_header(["Chic2ToKsKpiAlg/Chic2ToKsKpi.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_kskpi = Selection.new
sel_kskpi.select_track {
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     "==1"
            nChrn     "==1"
            nNet      "==0"
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=1"
         }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :pion, against: [:kaon, :proton]
            identify :kaon, against: [:pion, :proton]
         }
         .secondary_vertex_fit([:pip, :pim]) {
            build_virtual_particle(:K_S0).by_minimizing_mass_difference
            decay_length ">0.25"
            remove_used_particle_from_candidate_list
         }
         .kinematic_fit([:gamma, :K_S0, :kp, :pim]) {
            nominal
            constrain_four_momentum
            chi2_cut 60
         }

alg_kskpi.note(:ks_selection,
               "K_S0 candidates reconstructed from two oppositely charged tracks with loose vertex requirements (|Vz|<30 cm, Vr<10 cm) and no PID (both assumed pions); the candidate with mass closest to nominal M(K_S0) and secondary-vertex-fit decay length > 0.25 cm is chosen")
         .note(:chic2_signal_window,
               "|M(KS K pi) - M(chi_c2)| <= 15 MeV/c^2 signal region defined offline in ROOT")
         .with_decay_card(decay_card_kskpi)
         .apply(sel_kskpi)

alg_kskpi.execute_on([psip_data, psip_incMC, exMC_kskpi])

########################################################################
# Mode III: psi(3686) -> gamma pi+ pi- pi0
########################################################################
alg_3pi = Algorithm.new("Chic2ToPiPiPi0")
alg_3pi.set_header(["Chic2ToPiPiPi0Alg/Chic2ToPiPiPi0.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_3pi = Selection.new
sel_3pi.select_track {
          cos_theta 0.93
          Vz        10.0
          Vr        1.0
          nChrp     "==1"
          nChrn     "==1"
          nNet      "==0"
        }
       .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          nGam              ">=3"
       }
       .pid(method: :probability) {
          prob_cut 0.001
          identify :pion, against: [:kaon, :proton]
          npip ">=1"
          npim ">=1"
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 25
          npi0 ">=1"
       }
       .kinematic_fit([:gamma, :pip, :pim, :pi0]) {
          nominal
          constrain_four_momentum
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 60
          invariant_mass_of(:pip, :pim).within(0.0, 3.0)
       }

alg_3pi.note(:pi0_selection,
             "pi0 candidate reconstructed from two selected photons whose invariant mass is closest to pi0 nominal mass and satisfies |M_gg - M_pi0| < 10 MeV/c^2")
       .note(:pi0_recoil_mass_veto,
             "pi0 recoil mass < 3.0 GeV/c^2 required to suppress psi(3686) -> pi0 J/psi background")
       .note(:omega_veto,
             "M(gamma pi0) not in (0.7, 0.85) GeV/c^2 to veto psi(3686) -> omega pi+ pi- with omega -> gamma pi0")
       .note(:chic2_signal_window,
             "|M(pi+pi-pi0) - M(chi_c2)| <= 15 MeV/c^2 signal region defined offline in ROOT")
       .with_decay_card(decay_card_3pi)
       .apply(sel_3pi)

alg_3pi.execute_on([psip_data, psip_incMC, exMC_3pi])
