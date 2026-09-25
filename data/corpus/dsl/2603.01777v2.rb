# Search for psi(2S) -> Ds- pi+ + c.c. and psi(2S) -> Ds- rho+ + c.c.
# Two independent decay modes -> two Algorithm objects (Rule T1)

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # (2712.4 +/- 14.3)x10^6 psi(2S) events @ sqrt(s)=3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # Inclusive psi(2S) MC

# --- Decay card for signal Mode I: psi(2S) -> Ds- pi+ (Ds- -> phi e- nu_bar_e; phi -> K+ K-) ---
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.0000 D_s-  pi+                             VSS;
    Enddecay

    Decay D_s-
    1.0000 phi   e-   anti-nu_e                  PHOTOS ISGW2;
    Enddecay

    Decay phi
    1.0000 K+   K-                               VSS;
    Enddecay

    End
DECAYCARD

# --- Decay card for signal Mode II: psi(2S) -> Ds- rho+ (Ds- -> phi e- nu_bar_e; rho+ -> pi+ pi0; pi0 -> gamma gamma) ---
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 D_s-  rho+                            VVS_PWAVE 0.0 0.0 1.0 0.0 0.0 0.0;
    Enddecay

    Decay D_s-
    1.0000 phi   e-   anti-nu_e                  PHOTOS ISGW2;
    Enddecay

    Decay phi
    1.0000 K+   K-                               VSS;
    Enddecay

    Decay rho+
    1.0000 pi+  pi0                              VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                           PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC samples ---
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psi2S_to_Ds_pi_signal_mc"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end
exMC_modeI.save_to_config(format: :yaml, file_path: 'exMC_modeI_config')

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psi2S_to_Ds_rho_signal_mc"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end
exMC_modeII.save_to_config(format: :yaml, file_path: 'exMC_modeII_config')

######################################################################
### Mode I: psi(2S) -> Ds- pi+ ; Ds- -> phi e- nu_bar_e            ###
######################################################################
alg_modeI = Algorithm.new("Psi2S2DsPi")
alg_modeI.set_header(["Psi2S2DsPiAlg/Psi2S2DsPi.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_modeI = Selection.new
sel_modeI.select_track {
              cos_theta 0.93        # |cos(theta)| < 0.93
              Vr        1.0         # |Vxy| < 1 cm
              Vz        10.0        # |Vz| < 10 cm
              nChrp     ">=1"       # net charge = 0, K+ K- pi+ e- topology (2 positive, 2 negative)
              nChrn     ">=1"
              nNet      "==0"       # zero net charge
              nTot      "==4"       # exactly 4 charged tracks
            }
           .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14  # 700 ns TDC window
              angle_to_track    10.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
            }
           .pid(method: :probability) {
              prob_cut 0.0
              identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                              treat_as_electron_if_energy_above: 0.6
              identify :kaon, against: [:pion]           # K+/K- vs pi
              identify :pion, against: [:kaon, :electron] # pi+ vs K and e
              nkp  ">=1"
              nkm  ">=1"
              npip ">=1"
              nlm  ">=1"  # e- (via high-momentum leptons)
            }
           .partial_miss([5]) {
              # Miss the anti-nu_e in the decay chain; Ds- is inferred via recoil from pi+
              require_recoil_mass 1.91, 2.02   # Ds- recoil-mass window (~3 sigma), (1.91, 2.02) GeV/c^2
            }

alg_modeI
  .note(:phi_mass_window,
        "phi candidates selected by requiring M(K+K-) in (1.005, 1.035) GeV/c^2 (~3 sigma window around nominal phi mass)")
  .note(:electron_E_over_p,
        "electron candidates required to satisfy 0.86 < E/p < 1.03 to suppress e/pi misidentification")
  .note(:electron_pid_ratio,
        "electron PID: L(e) > 0.001 AND L(e)/(L(e)+L(pi)+L(K)) > 0.8")
  .note(:pmiss_cut,
        "|p_miss| > 0.02 GeV/c to suppress backgrounds without neutrinos; p_miss = p_psi(2S) - sum(p_visible)")
  .note(:Umiss_cut,
        "|U_miss| < 0.064 GeV; U_miss = E_miss - |p_miss| (neutrino identification)")
  .note(:Ds_recoil_method,
        "Ds- identified via recoil: E_Ds = E_psi(2S) - E_pi+, p_Ds = p_psi(2S) - p_pi+, M_Ds = sqrt(E_Ds^2 - |p_Ds|^2); signal window (1.91, 2.02) GeV/c^2")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)
alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])

######################################################################
### Mode II: psi(2S) -> Ds- rho+ ; rho+ -> pi+ pi0 ; pi0 -> gamma gamma
######################################################################
alg_modeII = Algorithm.new("Psi2S2DsRho")
alg_modeII.set_header(["Psi2S2DsRhoAlg/Psi2S2DsRho.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_modeII = Selection.new
sel_modeII.select_track {
              cos_theta 0.93
              Vr        1.0
              Vz        10.0
              nChrp     ">=1"
              nChrn     ">=1"
              nNet      "==0"
              nTot      "==4"
            }
           .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14
              angle_to_track    10.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              nGam              ">=2"   # need two photons for pi0
            }
           .pid(method: :probability) {
              prob_cut 0.0
              identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                              treat_as_electron_if_energy_above: 0.6
              identify :kaon, against: [:pion]
              identify :pion, against: [:kaon, :electron]
              nkp  ">=1"
              nkm  ">=1"
              npip ">=1"
              nlm  ">=1"
            }
           .kalman_kinematic_fit([:gamma, :gamma]) {
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              chi2_cut 200
              npi0 ">=1"
            }
           .partial_miss([5]) {
              # Miss anti-nu_e; Ds- inferred via recoil from rho+ (pi+ pi0)
              require_recoil_mass 1.89, 2.08   # Ds- recoil window for rho+ mode
            }

alg_modeII
  .note(:phi_mass_window,
        "phi candidates: M(K+K-) in (1.005, 1.035) GeV/c^2 (~3 sigma)")
  .note(:rho_mass_window,
        "rho+ candidates: M(pi+ pi0) in (0.61, 0.93) GeV/c^2 (~1 sigma)")
  .note(:electron_E_over_p,
        "electron: 0.86 < E/p < 1.03")
  .note(:electron_pid_ratio,
        "electron PID: L(e) > 0.001 AND L(e)/(L(e)+L(pi)+L(K)) > 0.8")
  .note(:pmiss_cut,
        "|p_miss| > 0.02 GeV/c")
  .note(:Umiss_cut,
        "|U_miss| < 0.10 GeV for Ds- rho+ mode")
  .note(:E_gamma_rest_cut,
        "E_gamma_rest < 0.12 GeV; sum energy of all photons NOT from the pi0 in rho+ decay; suppresses psi(2S) -> K+K-pi+pi-pi0 backgrounds")
  .note(:Ds_recoil_method,
        "Ds- via recoil against rho+: E_Ds = E_psi(2S) - E_rho+, M_Ds = sqrt(E_Ds^2 - |p_Ds|^2); signal window (1.89, 2.08) GeV/c^2")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)
alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])
