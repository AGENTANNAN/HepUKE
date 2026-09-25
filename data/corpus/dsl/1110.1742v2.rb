### Dataset description ###
psip_data   = DatasetManager.real_data.find("709_3686")     # psi(2S) real data (106 M events, 156.4 pb^-1)
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")  # psi(2S) inclusive MC

# Decay card for psi(2S) -> gamma chi_c2 -> gamma pi+ pi-
decay_card_pipi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2                            PHSP;
    Enddecay

    Decay chi_c2
    1.0000 pi+ pi-                                 PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for psi(2S) -> gamma chi_c2 -> gamma K+ K-
decay_card_kk = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2                            PHSP;
    Enddecay

    Decay chi_c2
    1.0000 K+ K-                                   PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples for the signal processes (phase space, used for normalization in PWA)
exMC_pipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_gamma_chic2_pipi"
  config.related_dataset = psip_data
  config.events         = 200000
  config.decay_card     = decay_card_pipi
  config.cross_section  = :default
end

exMC_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_gamma_chic2_kk"
  config.related_dataset = psip_data
  config.events         = 200000
  config.decay_card     = decay_card_kk
  config.cross_section  = :default
end

### Event selection (BOSS) — Mode I: psi(2S) -> gamma chi_c2 -> gamma pi+ pi- ###
alg_pipi = Algorithm.new("Chic2GamPiPi")
alg_pipi.set_header(["Chic2GamPiPiAlg/Chic2GamPiPi.h"])
        .set_constant({"ECMS" => [:double, 3.686]})

sel_pipi = Selection.new
sel_pipi.select_track {
           nChrp   "==1"    # exactly one positive track
           nChrn   "==1"    # exactly one negative track
           nNet    "==0"    # net charge zero
           cos_theta 0.93   # |cos(theta)| < 0.93
           Vz      10.0     # |Vz| < 10 cm
           Vr      1.0      # |Vr| < 1 cm
         }
         .select_photon {
           nGam                ">=1"    # at least one good photon (highest-E used as radiative gamma)
           energyThreshold_b   0.025    # >25 MeV in barrel (|cos theta|<0.8)
           energyThreshold_e   0.050    # >50 MeV in endcap (0.86<|cos theta|<0.92)
           angle_to_track      20.0     # >=20 deg away from any charged track
           tdc_emc_start       0
           tdc_emc_end         14
         }
         .pid(method: :probability) {
           prob_cut  0.001
           identify :pion, against: [:kaon, :proton]  # pi+/pi- PID: Prob(pi)>0.001
         }
         .kinematic_fit([:gamma, :pip, :pim]) {   # 4C fit under gamma pi+ pi- hypothesis
           nominal
           constrain_four_momentum
           chi2_cut 200                            # loose; tight cut (chi2_pi<60, chi2_pi<chi2_K) applied in ROOT
         }
         .kinematic_fit([:gamma, :kp, :km]) {     # competing gamma K+ K- hypothesis (stores chi2_K for ROOT-level pi/K arbitration)
           constrain_four_momentum
         }

alg_pipi.note(:electron_veto,
              "EMC deposited energy of each track < 1.4 GeV to remove e+e- -> (gamma)e+e- and psi(2S) -> (gamma)e+e-; " \
              "dE/dx within 3 sigma of expected value for each track; " \
              "for tracks in EMC-insensitive region 0.81<|cos theta|<0.86 the dE/dx window is tightened to 2 sigma.")
        .note(:muon_veto,
              "at least one charged track must have EMC deposited energy > 0.34 GeV to suppress " \
              "e+e- -> (gamma)mu+mu- and psi(2S) -> (gamma)mu+mu- (>99% of the mumu background removed).")
        .note(:barrel_endcap_gap,
              "photon showers in the barrel-endcap transition region (0.8<|cos theta|<0.86) are excluded.")

alg_pipi.with_decay_card(decay_card_pipi).apply(sel_pipi)
alg_pipi.execute_on([psip_data, psip_incMC, exMC_pipi])


### Event selection (BOSS) — Mode II: psi(2S) -> gamma chi_c2 -> gamma K+ K- ###
alg_kk = Algorithm.new("Chic2GamKK")
alg_kk.set_header(["Chic2GamKKAlg/Chic2GamKK.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

sel_kk = Selection.new
sel_kk.select_track {
         nChrp   "==1"
         nChrn   "==1"
         nNet    "==0"
         cos_theta 0.93
         Vz      10.0
         Vr      1.0
       }
       .select_photon {
         nGam                ">=1"
         energyThreshold_b   0.025
         energyThreshold_e   0.050
         angle_to_track      20.0
         tdc_emc_start       0
         tdc_emc_end         14
       }
       .pid(method: :probability) {
         prob_cut  0.001
         identify :kaon, against: [:pion, :proton]  # K+/K- PID: Prob(K)>Prob(pi) and Prob(K)>0.001
       }
       .kinematic_fit([:gamma, :kp, :km]) {   # 4C fit under gamma K+ K- hypothesis (nominal)
         nominal
         constrain_four_momentum
         chi2_cut 200                          # loose; tight arbitration (chi2_K<60, chi2_K<chi2_pi) done in ROOT
       }
       .kinematic_fit([:gamma, :pip, :pim]) { # competing gamma pi+ pi- hypothesis (stores chi2_pi for arbitration)
         constrain_four_momentum
       }

alg_kk.note(:electron_veto,
            "EMC deposited energy of each track < 1.4 GeV to remove e+e- -> (gamma)e+e- and psi(2S) -> (gamma)e+e-; " \
            "dE/dx within 3 sigma of expected value for each track.")
      .note(:barrel_endcap_gap,
            "photon showers in the barrel-endcap transition region (0.8<|cos theta|<0.86) are excluded.")

alg_kk.with_decay_card(decay_card_kk).apply(sel_kk)
alg_kk.execute_on([psip_data, psip_incMC, exMC_kk])
