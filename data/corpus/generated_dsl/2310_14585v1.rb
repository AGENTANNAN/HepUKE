# ==============================================================================
# psi(3686) -> Sigma+ Sigma_bar- omega / Sigma+ Sigma_bar- phi
#   Sigma+ -> p pi0,  Sigma_bar- -> p_bar pi0
#   omega  -> pi+ pi- pi0     (7C fit: 4C + 3 pi0 mass constraints)
#   phi    -> K+ K-           (2C fit: 4C + 1 pi0 mass constraint, 2nd pi0 missing)
# ==============================================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) real data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC

# --- Decay card: psi(3686) -> Sigma+ Sigma_bar- omega ---
decay_card_omega = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Sigma+ anti-Sigma- omega      PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0                        PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0                   PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0                   OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                   PHSP;
    Enddecay

    End
DECAYCARD

# --- Decay card: psi(3686) -> Sigma+ Sigma_bar- phi ---
decay_card_phi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Sigma+ anti-Sigma- phi        PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0                        PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0                   PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                         VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                   PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC for each decay mode (500k events) ---
exMC_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_sigmasigmabar_omega"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_omega
  config.cross_section   = :default
end

exMC_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_sigmasigmabar_phi"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_phi
  config.cross_section   = :default
end

# ==============================================================================
# Event selection: omega mode  psi(3686) -> Sigma+ Sigma_bar- omega
# ==============================================================================
alg_name_omega = "SigmaSigmaOmega"
alg_omega = Algorithm.new(alg_name_omega)
alg_omega.set_header(["#{alg_name_omega}Alg/#{alg_name_omega}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_omega = Selection.new
sel_omega.select_track {                     # charged track quality cuts
            cos_theta 0.93                   # |cos(theta)| < 0.93
            Vz        10.0                   # |Vz| < 10 cm
            Vr        2.0                    # Vr < 2 cm
            nChrp     ">=2"                  # at least 2 positive tracks
            nChrn     ">=2"                  # at least 2 negative tracks
          }
         .select_photon {                    # photon selection
            tdc_emc_start     0              # EMC timing 0-14
            tdc_emc_end       14
            angle_to_track    10.0           # min angle to any charged track (deg)
            energyThreshold_b 0.025          # barrel threshold 25 MeV
            energyThreshold_e 0.050          # endcap threshold 50 MeV
            nGam              ">=6"          # at least 6 photons (3 pi0)
          }
         .pid(method: :probability) {        # PID: probability method
            prob_cut 0.001
            identify :proton, against: [:pion, :kaon]   # p / p_bar
            identify :pion,   against: [:kaon]          # pi+ / pi-
            nprp ">=1"
            nprm ">=1"
            npip ">=1"
            npim ">=1"
          }
         .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C gamma gamma -> pi0
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 20
            npi0 ">=3"                       # at least three pi0
          }
         .kinematic_fit([:prp, :prm, :pip, :pim, :pi0, :pi0, :pi0]) {  # 7C fit
            nominal                          # 4C + three pi0 mass constraints
            constrain_four_momentum
            chi2_cut 45                      # loose cut; best combination kept
          }

# Post-fit Sigma mass window and background vetoes (handled on fit-corrected
# momenta in the ROOT analysis; captured here so they are not lost).
alg_omega.note(:sigma_mass_window,
               "Sigma+ mass window [1176, 1197] MeV/c^2 applied after the 7C fit on fit-corrected momenta")
         .note(:background_veto,
               "J/psi, eta and Lambda background vetoes applied after the 7C fit on fit-corrected momenta")

alg_omega.with_decay_card(decay_card_omega).apply(sel_omega)

# ==============================================================================
# Event selection: phi mode  psi(3686) -> Sigma+ Sigma_bar- phi
# ==============================================================================
alg_name_phi = "SigmaSigmaPhi"
alg_phi = Algorithm.new(alg_name_phi)
alg_phi.set_header(["#{alg_name_phi}Alg/#{alg_name_phi}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_phi = Selection.new
sel_phi.select_track {                       # charged track quality cuts
          cos_theta 0.93
          Vz        10.0
          Vr        2.0
          nChrp     ">=2"
          nChrn     ">=2"
        }
       .select_photon {                      # photon selection
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          nGam              ">=2"            # at least 2 photons
        }
       .pid(method: :probability) {          # PID: probability method
          prob_cut 0.001
          identify :proton, against: [:pion, :kaon]   # p / p_bar
          identify :kaon,   against: [:pion]          # K+ / K-
          nprp ">=1"
          nprm ">=1"
          nkp  ">=1"
          nkm  ">=1"
        }
       .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C gamma gamma -> pi0
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 20
          npi0 ">=1"                         # one reconstructed pi0 (best combo by mass)
        }
       .kinematic_fit([:prp, :prm, :kp, :km, :pi0]) {   # 2C fit
          nominal
          miss_track_of(:pi0)                # second pi0 from the other Sigma is missing
          constrain_four_momentum            # 4C + one pi0 mass constraint -> 2C
          chi2_cut 20
        }

# Post-fit Sigma mass window (handled on fit-corrected momenta in the ROOT analysis).
alg_phi.note(:sigma_mass_window,
             "Sigma+ mass window [1176, 1200] MeV/c^2 applied after the 2C fit on fit-corrected momenta")

alg_phi.with_decay_card(decay_card_phi).apply(sel_phi)

### Execute on datasets ###
root_files_omega = alg_omega.execute_on([psip_data, psip_incMC, exMC_omega])
root_files_phi   = alg_phi.execute_on([psip_data, psip_incMC, exMC_phi])