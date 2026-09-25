### Dataset description ###
# Total 2712.4e6 psi(3686) events collected in 2009, 2012, and 2021 at BESIII
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for the signal process: psi(3686) -> gamma chi_cJ -> gamma 2K+2K- omega (mode A)
decay_card_2K2Komega = <<~DECAYCARD
    Decay psi(3686)
    0.333  gamma  chi_c0                        PHOTOS  HELAMP 1.0 0.0 1.0 0.0;
    0.333  gamma  chi_c1                        PHOTOS  HELAMP 1.0 0.0 1.0 0.0 -1.0 0.0 -1.0 0.0;
    0.334  gamma  chi_c2                        PHOTOS  HELAMP 2.0 0.0 sqrt(1.5) 0.0 sqrt(0.5) 0.0 0.0 0.0 sqrt(0.5) 0.0 sqrt(1.5) 0.0 2.0 0.0;
    Enddecay

    Decay chi_c0
    1.000  K+ K+ K- K-  omega                   PHSP;
    Enddecay

    Decay chi_c1
    1.000  K+ K+ K- K-  omega                   PHSP;
    Enddecay

    Decay chi_c2
    1.000  K+ K+ K- K-  omega                   PHSP;
    Enddecay

    Decay omega
    1.000  pi+  pi-  pi0                        OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000  gamma gamma                          PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the second signal process: psi(3686) -> gamma chi_cJ -> gamma phi K+K- omega (mode B)
decay_card_phiKKomega = <<~DECAYCARD
    Alias  another_K+  K+
    Alias  another_K-  K-

    Decay psi(3686)
    0.333  gamma  chi_c0                        PHOTOS  HELAMP 1.0 0.0 1.0 0.0;
    0.333  gamma  chi_c1                        PHOTOS  HELAMP 1.0 0.0 1.0 0.0 -1.0 0.0 -1.0 0.0;
    0.334  gamma  chi_c2                        PHOTOS  HELAMP 2.0 0.0 sqrt(1.5) 0.0 sqrt(0.5) 0.0 0.0 0.0 sqrt(0.5) 0.0 sqrt(1.5) 0.0 2.0 0.0;
    Enddecay

    Decay chi_c0
    1.000  phi  another_K+  another_K-  omega   PHSP;
    Enddecay

    Decay chi_c1
    1.000  phi  another_K+  another_K-  omega   PHSP;
    Enddecay

    Decay chi_c2
    1.000  phi  another_K+  another_K-  omega   PHSP;
    Enddecay

    Decay phi
    1.000  K+  K-                               VSS;
    Enddecay

    Decay omega
    1.000  pi+  pi-  pi0                        OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000  gamma gamma                          PHSP;
    Enddecay

    End
DECAYCARD

exMC_2K2Kw = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "chi_cJ_to_2K2K_omega"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_2K2Komega
  config.cross_section   = :default
end
exMC_2K2Kw.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_phiKKw = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "chi_cJ_to_phiKK_omega"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_phiKKomega
  config.cross_section   = :default
end
exMC_phiKKw.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
# Both signal processes share the same final state gamma 2K+ 2K- pi+ pi- pi0 (6 charged tracks,
# net charge 0, at least 3 photons); a single Algorithm covers both.
alg_name = "chi_cJ_2K2Komega"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.686] })

event_selection = Selection.new
event_selection.select_track {                       # Charged track selection
                 cos_theta 0.93                       # |cos(theta)| < 0.93
                 Vz        10.0                       # |Vz| < 10 cm along beam axis
                 Vr        1.0                        # |Vxy| < 1 cm in transverse plane
                 nChrp     "==3"                      # 3 positive charged tracks (2K+ + pi+)
                 nChrn     "==3"                      # 3 negative charged tracks (2K- + pi-)
                 nNet      "==0"                      # net charge = 0
               }
               .select_photon {                       # Photon selection
                 tdc_emc_start   0
                 tdc_emc_end     700
                 angle_to_track  10.0                 # min angle to nearest charged track (deg)
                 energyThreshold_b 0.025              # E > 25 MeV in EMC barrel
                 energyThreshold_e 0.050              # E > 50 MeV in EMC end-cap
                 nGam            ">=3"                # at least 3 photons (1 transition gamma + 2 for pi0)
               }
               .pid(method: :probability) {          # Kaon / pion PID via dE/dx + TOF
                 prob_cut 0.001
                 identify :kaon, against: [:pion]     # L(K) > L(pi)  -> kaon
                 identify :pion, against: [:kaon]     # L(pi) > L(K)  -> pion
                 nkp ">=2"; nkm ">=2"                 # at least 2 K+ and 2 K-
                 npip ">=1"; npim ">=1"               # at least 1 pi+ and 1 pi-
               }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct pi0 from photon pair (1-C)
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 chi2_cut 25
                 npi0 ">=1"
               }
               # 5C kinematic fit: 4-momentum constraint on psi(3686) + pi0 mass constraint
               # (5C in the paper). We take :pi0 as a mass-constrained participant, so the
               # remaining fit is 4C on the total four-momentum.
               .kinematic_fit([:gamma, :kp, :kp, :km, :km, :pip, :pim, :pi0]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 40                          # events with chi2_5C < 40 retained
               }

alg.with_decay_card(decay_card_2K2Komega).apply(event_selection)
root_files = alg.execute_on([psip_data, psip_incMC, exMC_2K2Kw, exMC_phiKKw])
