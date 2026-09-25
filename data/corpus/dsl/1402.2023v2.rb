# Analysis: chi_cJ -> eta' K+ K- (J=1,2) from psi(3686) -> gamma chi_cJ
# Two eta' decay modes: I) eta' -> gamma rho^0 (rho^0 -> pi+pi-); II) eta' -> eta pi+ pi-, eta -> gamma gamma

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for mode I: eta' -> gamma rho0
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma chi_c1                       PHOTOS HELAMP 1.0 0.0 1.0 0.0 -1.0 0.0 -1.0 0.0;
    Enddecay

    Decay chi_c1
    1.0000  eta' K+ K-                         PHSP;
    Enddecay

    Decay eta'
    1.0000  gamma rho0                         PHSP;
    Enddecay

    Decay rho0
    1.0000  pi+ pi-                            VSS;
    Enddecay

    End
DECAYCARD

# Decay card for mode II: eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma chi_c1                       PHOTOS HELAMP 1.0 0.0 1.0 0.0 -1.0 0.0 -1.0 0.0;
    Enddecay

    Decay chi_c1
    1.0000  eta' K+ K-                         PHSP;
    Enddecay

    Decay eta'
    1.0000  eta pi+ pi-                        PHSP;
    Enddecay

    Decay eta
    1.0000  gamma gamma                        PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "chic1_etap_KK_modeI"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "chic1_etap_KK_modeII"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Mode I algorithm: eta' -> gamma rho0 ###
alg_modeI = Algorithm.new("ChiCJEtapKKModeI")
alg_modeI.set_header(["ChiCJEtapKKModeIAlg/ChiCJEtapKKModeI.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_modeI = Selection.new
sel_modeI.select_track {
             cos_theta 0.93          # |cos(theta)| < 0.93
             Vz        10.0          # |Vz| < 10 cm
             Vr        1.0           # |Vxy| < 1 cm
             nChrp     "==2"         # 2 positive tracks (K+ and pi+)
             nChrn     "==2"         # 2 negative tracks (K- and pi-)
             nNet      "==0"         # net charge zero
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             angle_to_track    5.0    # >5 degrees from any charged track
             energyThreshold_b 0.025  # barrel: E > 25 MeV
             energyThreshold_e 0.050  # endcap: E > 50 MeV
             nGam              ">=2"  # at least 2 good photons for mode I
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]   # Prob(K)>Prob(pi), Prob(K)>Prob(p)
             nkp ">=1"
             nkm ">=1"
           }
          .assign({:chrgp => :pip, :chrgn => :pim})       # remaining tracks assumed to be pions
          .remove([:kp <= :pip, :km <= :pim])
          .kinematic_fit([:gamma, :gamma, :kp, :km, :pip, :pim]) {
             nominal
             constrain_four_momentum
             chi2_cut 40   # chi2_4C < 40 (mode I)
           }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.note(:eta_prime_photon_assignment,
               "photon with smaller |M(gamma pi+ pi-) - M(eta')| assigned as the eta' daughter photon; other one tagged as the radiative photon from psi(3686)")
         .note(:pi0_veto,        "|M(gamma gamma) - M(pi0)| > 15 MeV/c^2 to remove pi0 background")
         .note(:jpsi_veto_pipi,  "|M(pi+ pi-)_rec - M(J/psi)| > 8 MeV/c^2 to reject psi(3686) -> pi+ pi- J/psi backgrounds")
         .note(:jpsi_veto_gg,    "|M(gamma gamma)_rec - M(J/psi)| > 22 MeV/c^2 to reject psi(3686) -> gamma gamma J/psi backgrounds")
         .note(:etap_mass_window,"|M(gamma pi+ pi-) - M(eta')| < 15 MeV/c^2 for eta' signal selection")
         .note(:chicJ_mass_window,"chi_cJ mass windows for J=0/1/2: 30/15/16 MeV/c^2 around M(chi_cJ) on M(gamma pi+ pi- K+ K-)")

alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])

### Mode II algorithm: eta' -> eta pi+ pi-, eta -> gamma gamma ###
alg_modeII = Algorithm.new("ChiCJEtapKKModeII")
alg_modeII.set_header(["ChiCJEtapKKModeIIAlg/ChiCJEtapKKModeII.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_modeII = Selection.new
sel_modeII.select_track {
              cos_theta 0.93
              Vz        10.0
              Vr        1.0
              nChrp     "==2"
              nChrn     "==2"
              nNet      "==0"
            }
           .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14
              angle_to_track    5.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              nGam              ">=3"   # at least 3 good photons for mode II
            }
           .pid(method: :probability) {
              prob_cut 0.001
              identify :kaon, against: [:pion, :proton]
              nkp ">=1"
              nkm ">=1"
            }
           .assign({:chrgp => :pip, :chrgn => :pim})
           .remove([:kp <= :pip, :km <= :pim])
           .kinematic_fit([:gamma, :gamma, :gamma, :kp, :km, :pip, :pim]) {
              nominal
              constrain_four_momentum
              chi2_cut 50   # chi2_4C < 50 (mode II)
            }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.note(:pi0_veto,          "reject events if any gamma-gamma pair has |M(gamma gamma) - M(pi0)| < 20 MeV/c^2")
          .note(:eta_selection,     "eta candidate = photon pair with M(gamma gamma) closest to M(eta); |M(gamma gamma) - M(eta)| < 25 MeV/c^2")
          .note(:etap_mass_window,  "|M(eta pi+ pi-) - M(eta')| < 25 MeV/c^2 for eta' signal")
          .note(:chicJ_mass_window, "chi_cJ mass windows for J=0/1/2: 36/18/18 MeV/c^2 on M(gamma gamma pi+ pi- K+ K-)")

alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])
