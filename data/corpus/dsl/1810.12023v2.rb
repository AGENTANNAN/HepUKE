### dataset description ###
psi3686_data = DatasetManager.real_data.find("709_3686")
psi3686_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Common: all modes use psi(3686) -> pi0 h_c, h_c -> hadrons
# The tag pi0 is the lower-energy pi0 in events with multiple pi0s

### Mode I: h_c -> p pbar pi+ pi- ###
decay_card_modeI = <<~DECAYCARD
    Decay psi(3686)
    1.0000 pi0 h_c    PHSP;
    Enddecay

    Decay h_c
    1.0000 p+ anti-p- pi+ pi-    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_psip_pi0_hc_pppipi"
    config.related_dataset = psi3686_data
    config.events = 100000
    config.decay_card = decay_card_modeI
    config.cross_section = :default
end

alg_modeI = Algorithm.new("HcToPpbarPiPi")
alg_modeI.set_header(["HcToPpbarPiPiAlg/HcToPpbarPiPi.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_modeI = Selection.new
sel_modeI.select_track {
              cos_theta  0.93
              Vz         10.0
              Vr         1.0
              nChrp      "==2"
              nChrn      "==2"
              nNet       "==0"
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
              identify :proton, against: [:pion, :kaon]
              nprp ">=1"; nprm ">=1"
          }
          .remove([:prp <= :chrgp, :prm <= :chrgn])
          .assign({chrgp: :pip, chrgn: :pim})
          .kalman_kinematic_fit([:gamma, :gamma]) {
              build_virtual_particle(:pi0).by_minimizing_mass_difference
          }
          .kinematic_fit([:pi0, :prp, :prm, :pip, :pim]) {
              nominal
              vertex_fit([1, 2, 3, 4])
              constrain_four_momentum
              chi2_cut 200
          }

alg_modeI
  .note(:tag_pi0, "tag pi0 from psi(3686)->pi0 h_c identified via recoil mass; h_c signal region |RM(pi0)-3.525| < 8 MeV/c^2; for events with multiple pi0 candidates, use lower-energy pi0 as tag")
  .note(:proton_pid, "proton/antiproton identified via TOF+dE/dx probability; proton photon angle cut >20 deg applied in ROOT (DSL uses looser 10 deg for all tracks)")
  .note(:competing_hypothesis_veto, "chi2_4C_exp(2gamma) < chi2_4C_unexp(1gamma) to suppress psi(3686)->gamma chi_c2 backgrounds; expectation from 4-momentum fit with expected vs unexpected photon count")
  .note(:excess_photon_combinations, "if >2 photon candidates, loop all gamma-gamma combinations and keep best chi2_4C; best pi0 chosen by minimizing |M(gamma gamma)-M_pi0|")
  .note(:mass_windows, "post-fit mass windows (Table I): |RM(pi+pi-)-M_J/psi|>18, |M(pi+pi-pi0)-M_eta|>14, |M(pi+pi-pi0)-M_omega|>6 MeV/c^2 to suppress psi(3686)->pi+pi-J/psi, pi0 eta, pi0 omega backgrounds")
  .note(:signal_extraction, "h_c yield from unbinned ML fit to RM(pi0) spectrum; signal = MC shape convoluted with Gaussian (data-MC resolution); background = ARGUS function; h_c->ppbar pi+pi- observed at 7.4 sigma; BF = (2.89+/-0.32+/-0.55)x10^{-3}")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

### Mode II: h_c -> pi+ pi- pi0 ###
decay_card_modeII = <<~DECAYCARD
    Decay psi(3686)
    1.0000 pi0 h_c    PHSP;
    Enddecay

    Decay h_c
    1.0000 pi+ pi- pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_psip_pi0_hc_3pi"
    config.related_dataset = psi3686_data
    config.events = 100000
    config.decay_card = decay_card_modeII
    config.cross_section = :default
end

alg_modeII = Algorithm.new("HcToPiPiPi0")
alg_modeII.set_header(["HcToPiPiPi0Alg/HcToPiPiPi0.h"])
           .set_constant({"ECMS" => [:double, 3.686]})

sel_modeII = Selection.new
sel_modeII.select_track {
               cos_theta  0.93
               Vz         10.0
               Vr         1.0
               nChrp      "==1"
               nChrn      "==1"
               nNet       "==0"
           }
           .select_photon {
               tdc_emc_start     0
               tdc_emc_end       14
               energyThreshold_b 0.025
               energyThreshold_e 0.050
               angle_to_track    10.0
               nGam              ">=4"
           }
           .assign({chrgp: :pip, chrgn: :pim})
           .kalman_kinematic_fit([:gamma, :gamma]) {
               build_virtual_particle(:pi0_tag).by_minimizing_mass_difference
           }
           .kalman_kinematic_fit([:gamma, :gamma]) {
               build_virtual_particle(:pi0_hc).by_minimizing_mass_difference
           }
           .kinematic_fit([:pi0_tag, :pi0_hc, :pip, :pim]) {
               nominal
               vertex_fit([2, 3])
               constrain_four_momentum
               chi2_cut 200
           }

alg_modeII
  .note(:tag_pi0, "tag pi0 from psi(3686)->pi0 h_c (lower-energy pi0); h_c signal region |RM(pi0_tag)-3.525| < 8 MeV/c^2")
  .note(:competing_hypothesis_veto, "chi2_4C_exp(4gamma) < chi2_4C_unexp(3gamma) to suppress backgrounds with different photon multiplicity")
  .note(:excess_photon_combinations, "if >4 photon candidates, loop combinations and keep best chi2_4C; pi0 mass constraint 107<M(gamma gamma)<163 MeV/c^2 applied in kinematic fit")
  .note(:mass_windows, "post-fit mass windows (Table I): |RM(pi+pi-)-M_J/psi|>74, |RM(pi0)-M_omega|>32 MeV/c^2 to suppress psi(3686)->pi+pi-J/psi and pi0 omega backgrounds")
  .note(:signal_extraction, "h_c yield from unbinned ML fit to RM(pi0_tag); h_c->pi+pi-pi0 observed at 4.6 sigma; BF = (1.60+/-0.40+/-0.32)x10^{-3}")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

### Mode III: h_c -> 2(pi+ pi-) pi0 ###
decay_card_modeIII = <<~DECAYCARD
    Decay psi(3686)
    1.0000 pi0 h_c    PHSP;
    Enddecay

    Decay h_c
    1.0000 pi+ pi- pi+ pi- pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_psip_pi0_hc_4pipi0"
    config.related_dataset = psi3686_data
    config.events = 100000
    config.decay_card = decay_card_modeIII
    config.cross_section = :default
end

alg_modeIII = Algorithm.new("HcTo4PiPi0")
alg_modeIII.set_header(["HcTo4PiPi0Alg/HcTo4PiPi0.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

sel_modeIII = Selection.new
sel_modeIII.select_track {
                cos_theta  0.93
                Vz         10.0
                Vr         1.0
                nChrp      "==2"
                nChrn      "==2"
                nNet       "==0"
            }
            .select_photon {
                tdc_emc_start     0
                tdc_emc_end       14
                energyThreshold_b 0.025
                energyThreshold_e 0.050
                angle_to_track    10.0
                nGam              ">=4"
            }
            .assign({chrgp: :pip, chrgn: :pim})
            .kalman_kinematic_fit([:gamma, :gamma]) {
                build_virtual_particle(:pi0_tag).by_minimizing_mass_difference
            }
            .kalman_kinematic_fit([:gamma, :gamma]) {
                build_virtual_particle(:pi0_hc).by_minimizing_mass_difference
            }
            .kinematic_fit([:pi0_tag, :pi0_hc, :pip, :pim, :pip, :pim]) {
                nominal
                vertex_fit([2, 3, 4, 5])
                constrain_four_momentum
                chi2_cut 200
            }

alg_modeIII
  .note(:tag_pi0, "tag pi0 from psi(3686)->pi0 h_c (lower-energy pi0); h_c signal region |RM(pi0_tag)-3.525| < 8 MeV/c^2")
  .note(:competing_hypothesis_veto, "chi2_4C_exp(4gamma) < chi2_4C_unexp(3gamma) to suppress backgrounds with different photon multiplicity")
  .note(:excess_photon_combinations, "if >4 photon candidates, loop combinations and keep best chi2_4C")
  .note(:mass_windows, "post-fit mass windows (Table I): |RM(pi+pi-)-M_J/psi|>20, |RM(pi+pi-)-M_J/psi|>22, |M(pi+pi-pi0)-M_eta|>16, |M(pi+pi-pi0)-M_omega|>20 MeV/c^2")
  .note(:signal_extraction, "h_c yield from unbinned ML fit to RM(pi0_tag); h_c->2(pi+pi-)pi0 observed at 9.1 sigma; BF = (7.44+/-0.94+/-1.52)x10^{-3}")
  .note(:intermediate_resonances, "rho0 peaks observed in pi+pi- invariant mass projections; systematic uncertainties from physics model evaluated with alternative simulation samples incorporating intermediate states")
  .with_decay_card(decay_card_modeIII)
  .apply(sel_modeIII)

### Mode IV: h_c -> 3(pi+ pi-) pi0 (upper limit) ###
decay_card_modeIV = <<~DECAYCARD
    Decay psi(3686)
    1.0000 pi0 h_c    PHSP;
    Enddecay

    Decay h_c
    1.0000 pi+ pi- pi+ pi- pi+ pi- pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeIV = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_psip_pi0_hc_6pipi0"
    config.related_dataset = psi3686_data
    config.events = 100000
    config.decay_card = decay_card_modeIV
    config.cross_section = :default
end

alg_modeIV = Algorithm.new("HcTo6PiPi0")
alg_modeIV.set_header(["HcTo6PiPi0Alg/HcTo6PiPi0.h"])
           .set_constant({"ECMS" => [:double, 3.686]})

sel_modeIV = Selection.new
sel_modeIV.select_track {
               cos_theta  0.93
               Vz         10.0
               Vr         1.0
               nChrp      "==3"
               nChrn      "==3"
               nNet       "==0"
           }
           .select_photon {
               tdc_emc_start     0
               tdc_emc_end       14
               energyThreshold_b 0.025
               energyThreshold_e 0.050
               angle_to_track    10.0
               nGam              ">=4"
           }
           .assign({chrgp: :pip, chrgn: :pim})
           .kalman_kinematic_fit([:gamma, :gamma]) {
               build_virtual_particle(:pi0_tag).by_minimizing_mass_difference
           }
           .kalman_kinematic_fit([:gamma, :gamma]) {
               build_virtual_particle(:pi0_hc).by_minimizing_mass_difference
           }
           .kinematic_fit([:pi0_tag, :pi0_hc, :pip, :pim, :pip, :pim, :pip, :pim]) {
               nominal
               vertex_fit([2, 3, 4, 5, 6, 7])
               constrain_four_momentum
               chi2_cut 200
           }

alg_modeIV
  .note(:tag_pi0, "tag pi0 from psi(3686)->pi0 h_c (lower-energy pi0); h_c signal region |RM(pi0_tag)-3.525| < 8 MeV/c^2")
  .note(:competing_hypothesis_veto, "chi2_4C_exp(4gamma) < chi2_4C_unexp(3gamma)")
  .note(:excess_photon_combinations, "if >4 photon candidates, loop combinations and keep best chi2_4C")
  .note(:mass_windows, "post-fit mass windows (Table I): |RM(pi+pi- pi0_hc lowE)-M_J/psi|>18, |RM(pi+pi-)-M_J/psi|>20, |M(pi+pi-pi0)-M_eta|>16 MeV/c^2")
  .note(:signal_extraction, "no significant signal; 90% CL Bayesian upper limit; BF < 8.7x10^{-3} at 90% CL; largest upper limit taken from different fitting models/ranges")
  .with_decay_card(decay_card_modeIV)
  .apply(sel_modeIV)

### Mode V: h_c -> K+ K- pi+ pi- (upper limit) ###
decay_card_modeV = <<~DECAYCARD
    Decay psi(3686)
    1.0000 pi0 h_c    PHSP;
    Enddecay

    Decay h_c
    1.0000 K+ K- pi+ pi-    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeV = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_psip_pi0_hc_KKpipi"
    config.related_dataset = psi3686_data
    config.events = 100000
    config.decay_card = decay_card_modeV
    config.cross_section = :default
end

alg_modeV = Algorithm.new("HcToKKPiPi")
alg_modeV.set_header(["HcToKKPiPiAlg/HcToKKPiPi.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_modeV = Selection.new
sel_modeV.select_track {
              cos_theta  0.93
              Vz         10.0
              Vr         1.0
              nChrp      "==2"
              nChrn      "==2"
              nNet       "==0"
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
              identify :kaon, against: [:pion, :proton]
              nkp ">=1"; nkm ">=1"
          }
          .remove([:kp <= :chrgp, :km <= :chrgn])
          .assign({chrgp: :pip, chrgn: :pim})
          .kalman_kinematic_fit([:gamma, :gamma]) {
              build_virtual_particle(:pi0).by_minimizing_mass_difference
          }
          .kinematic_fit([:pi0, :kp, :km, :pip, :pim]) {
              nominal
              vertex_fit([1, 2, 3, 4])
              constrain_four_momentum
              chi2_cut 200
          }

alg_modeV
  .note(:tag_pi0, "tag pi0 from psi(3686)->pi0 h_c; h_c signal region |RM(pi0)-3.525| < 8 MeV/c^2")
  .note(:kaon_pid, "kaon identified via TOF+dE/dx probability; remaining tracks assigned as pions")
  .note(:competing_hypothesis_veto, "chi2_4C_exp(2gamma) < chi2_4C_unexp(1gamma) to suppress psi(3686)->gamma chi_c2 backgrounds; chi_c2->K+K-pi+pi- peaking background included as additional component in fit")
  .note(:excess_photon_combinations, "if >2 photon candidates, loop combinations and keep best chi2_4C")
  .note(:mass_windows, "post-fit mass windows (Table I): |RM(pi+pi-)-M_J/psi|>22, |M(pi+pi-pi0)-M_eta|>16, |M(pi+pi-pi0)-M_omega|>20 MeV/c^2")
  .note(:signal_extraction, "no significant signal; 90% CL Bayesian upper limit; BF < 5.8x10^{-4} at 90% CL; additional chi_c2->K+K-pi+pi- background component included in fit")
  .note(:chi_c2_background, "MC study shows peaking background from psi(3686)->gamma chi_c2, chi_c2->K+K-pi+pi- contributes to this mode; included as additional component in RM(pi0) fit")
  .with_decay_card(decay_card_modeV)
  .apply(sel_modeV)

alg_modeI.execute_on([psi3686_data, psi3686_incMC, exMC_modeI])
alg_modeII.execute_on([psi3686_data, psi3686_incMC, exMC_modeII])
alg_modeIII.execute_on([psi3686_data, psi3686_incMC, exMC_modeIII])
alg_modeIV.execute_on([psi3686_data, psi3686_incMC, exMC_modeIV])
alg_modeV.execute_on([psi3686_data, psi3686_incMC, exMC_modeV])