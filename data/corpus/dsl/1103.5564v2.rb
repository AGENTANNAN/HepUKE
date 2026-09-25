# BESIII paper 1103.5564v2 — chi_cJ -> gamma V (V = phi, rho0, omega)
# psi' -> gamma_l chi_cJ -> gamma_l gamma_h V, with V -> K+K-, pi+pi-, pi+pi-pi0.
# Samples: (1.06 +- 0.04) x 10^8 psi' events + ~42.6 pb^-1 continuum data at 3.65 GeV.
# The three vector-meson channels have different final states and different
# kinematic-fit hypotheses, so each gets its own Algorithm + Selection chain.

### Dataset preparation ###
psip_data   = DatasetManager.real_data.find("709_3686")      # psi(2S) data
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")   # inclusive psi' MC
cont_data   = DatasetManager.real_data.find("709_3650")      # continuum data at 3.65 GeV
cont_incMC  = DatasetManager.inclusive_mc.find("709_3650")   # continuum inclusive MC

# psi' -> gamma gamma phi, phi -> K+ K-
decay_card_gammaphi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma gamma phi          PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                    VSS;
    Enddecay

    End
DECAYCARD

# psi' -> gamma gamma rho0, rho0 -> pi+ pi-
decay_card_gammarho = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma gamma rho0         PHSP;
    Enddecay

    Decay rho0
    1.0000 pi+ pi-                  VSS;
    Enddecay

    End
DECAYCARD

# psi' -> gamma gamma omega, omega -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_gammaomega = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma gamma omega        PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0              OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma              PHSP;
    Enddecay

    End
DECAYCARD

exMC_gammaphi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_gamma_phi_KK"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_gammaphi
    config.cross_section   = :default
end

exMC_gammarho = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_gamma_rho0_pipi"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_gammarho
    config.cross_section   = :default
end

exMC_gammaomega = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_gamma_omega_pipipi0"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_gammaomega
    config.cross_section   = :default
end

### Event selection (BOSS) — Mode 1: psi' -> gamma gamma phi, phi -> K+ K- ###
alg_name_phi = "ChicJGammaPhi"
alg_phi = Algorithm.new(alg_name_phi)
alg_phi.set_header(["#{alg_name_phi}Alg/#{alg_name_phi}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_phi = Selection.new
sel_phi.select_track {
          cos_theta 0.93   # |cos(theta)| < 0.93
          Vz        10.0   # within +-10 cm of the IP along the beam direction
          Vr         1.0   # within 1 cm of the beam line in the transverse plane
          nChrp    "==1"   # exactly two charged tracks with net charge zero
          nChrn    "==1"
          nNet     "==0"
        }
        .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0    # photon separated by >= 10 deg from any charged track
          energyThreshold_b 0.025   # barrel (|cos theta| < 0.8) : E > 25 MeV
          energyThreshold_e 0.050   # end-cap (0.86 < |cos theta| < 0.92) : E > 50 MeV
          nGam              ">=2"   # two radiative photons
        }
        .pid(method: :probability) {
          prob_cut 0.001
          # Only the lower-momentum charged track is required to be identified,
          # because PID efficiency is lower for p > 1 GeV/c tracks.
          identify :kaon, against: [:pion, :proton]
          nkp ">=1"
          nkm ">=1"
        }
        # 4C kinematic fit under energy-momentum conservation; the candidate with the
        # smallest chi2 is kept when the event has more than one combination.
        .kinematic_fit([:gamma, :gamma, :kp, :km]) {
          nominal
          constrain_four_momentum
          chi2_cut 200   # loose BOSS cut; the paper requires chi2_4C <= 100 (applied in ROOT)
        }

alg_phi.note(:phi_signal_and_sideband,
             "phi signal region |M(K+K-) - M_phi| <= 0.01 GeV/c^2; sideband " \
             "1.05 <= M(K+K-) <= 1.07 GeV/c^2. Windows applied in ROOT.")
       .note(:background_veto,
             "Background from multi-photon hadronic psi' decays suppressed by requiring " \
             "|M(gamma_l gamma_h) - M_eta| >= 25 MeV/c^2. Also the psi' -> gamma_h eta', " \
             "eta' -> gamma_l V background is suppressed by |M(gamma_l V) - M_eta'| > 15 MeV/c^2. " \
             "Applied in ROOT after the 4C fit.")
       .note(:mass_spectrum_fit,
             "The gamma_h phi invariant mass distribution is fitted with MC signal shapes for the " \
             "three chi_cJ resonances plus a background composed of the phi mass sideband " \
             "distribution and a 2nd-order polynomial. Fitted yields: " \
             "15.0 +- 6.6 (chi_c0, <16.2e-6 at 90% C.L.), 42.6 +- 8.6 (chi_c1, 6.4 sigma), " \
             "4.6 +- 4.9 (chi_c2, <8.1e-6). Efficiencies 32.4% / 34.6% / 32.6%.")
       .note(:polarization_fit,
             "Transverse polarization fraction f_T extracted from a likelihood fit to the " \
             "cos(Theta) distribution in the chi_c1 signal region " \
             "3.49 <= M(gamma_h V) <= 3.52 GeV/c^2, where Theta is the angle between the vector " \
             "meson flight direction in the chi_c1 rest frame and the K+ direction in the phi " \
             "rest frame. f_T = 0.29 +0.13+0.10 -0.12-0.09 for chi_c1 -> gamma phi. ROOT-level.")
       .note(:systematics,
             "Tracking 4.0%, PID 2.0%, photon detection 2.0%, 4C fit 0.7%, selection-efficiency " \
             "uncertainty 2.0% (eta veto 1.9%, phi selection 0.5%), background shape, binning, " \
             "fit range and sideband regions, signal shape, N_psi' 3.8% and " \
             "B(psi' -> gamma chi_cJ) / B(phi -> K+K-) from the PDG. Total 8.8-9.3%.")

alg_phi.with_decay_card(decay_card_gammaphi).apply(sel_phi)
root_files_phi = alg_phi.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_gammaphi])

### Event selection (BOSS) — Mode 2: psi' -> gamma gamma rho0, rho0 -> pi+ pi- ###
alg_name_rho = "ChicJGammaRho0"
alg_rho = Algorithm.new(alg_name_rho)
alg_rho.set_header(["#{alg_name_rho}Alg/#{alg_name_rho}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_rho = Selection.new
sel_rho.select_track {
          cos_theta 0.93
          Vz        10.0
          Vr         1.0
          nChrp    "==1"
          nChrn    "==1"
          nNet     "==0"
        }
        .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          nGam              ">=2"
        }
        .pid(method: :probability) {
          prob_cut 0.001
          identify :pion, against: [:kaon, :proton]   # only the lower-momentum track must be identified
          npip ">=1"
          npim ">=1"
        }
        .kinematic_fit([:gamma, :gamma, :pip, :pim]) {
          nominal
          constrain_four_momentum
          chi2_cut 200   # loose BOSS cut; the paper requires chi2_4C <= 100 (applied in ROOT)
        }

alg_rho.note(:rho_signal_and_sideband,
             "rho0 signal region |M(pi+pi-) - M_rho| <= 0.2 GeV/c^2; sideband " \
             "1.25 <= M(pi+pi-) <= 1.65 GeV/c^2. Windows applied in ROOT.")
       .note(:background_veto,
             "QED backgrounds e+e- -> gamma e+e- and gamma mu+mu- with the leptons misidentified as " \
             "pions are rejected by requiring EMC energy over MDC momentum E_EMC/(c p_MDC) < 0.8 " \
             "for each track (electron rejection), and by removing tracks with more than three " \
             "muon-chamber layers with hits (muon rejection). The remaining QCD background is " \
             "removed by requiring the pion opening angle cos(theta_pi+pi-) > -0.8 and the photon " \
             "opening angle -0.98 < cos(theta_gamma_l gamma_h) < 0.5 in the laboratory frame. " \
             "In addition |M(gamma_l gamma_h) - M_eta| >= 25 MeV/c^2 and " \
             "M(gamma_l gamma_h) >= 600 MeV/c^2 are required. Not expressible in the DSL " \
             "selection blocks; applied in ROOT.")
       .note(:mass_spectrum_fit,
             "The gamma_h rho0 invariant mass distribution is fitted with MC signal shapes for the " \
             "three chi_cJ resonances plus a background composed of the rho0 mass sideband " \
             "distribution and a 2nd-order polynomial. Fitted yields: 6 +- 12 (chi_c0, <10.5e-6), " \
             "432 +- 25 (chi_c1, >10 sigma), 13 +- 11 (chi_c2, <20.8e-6). " \
             "Efficiencies 22.6% / 19.4% / 15.7%.")
       .note(:polarization_fit,
             "Transverse polarization fraction f_T from a likelihood fit to the cos(Theta) " \
             "distribution, Theta being the angle between the rho0 flight direction in the chi_c1 " \
             "rest frame and the pi+ direction in the rho0 rest frame; " \
             "f_T = 0.158 +- 0.034 +0.015 -0.014 for chi_c1 -> gamma rho0. ROOT-level.")
       .note(:systematics,
             "Tracking 4.0%, PID 2.0%, photon detection 2.0%, 4C fit 0.7%, selection-efficiency " \
             "uncertainty 5.0% (electron/muon rejection 4.5%, cos(theta_gamma gamma) 0.9%, " \
             "pi0/eta veto 1.2%, rho0 selection 1.4%), background shape, binning, fit range and " \
             "sideband regions, signal shape, N_psi' and PDG branching fractions. Total ~8.3%.")

alg_rho.with_decay_card(decay_card_gammarho).apply(sel_rho)
root_files_rho = alg_rho.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_gammarho])

### Event selection (BOSS) — Mode 3: psi' -> gamma gamma omega, omega -> pi+ pi- pi0 ###
alg_name_omega = "ChicJGammaOmega"
alg_omega = Algorithm.new(alg_name_omega)
alg_omega.set_header(["#{alg_name_omega}Alg/#{alg_name_omega}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_omega = Selection.new
sel_omega.select_track {
            cos_theta 0.93
            Vz        10.0
            Vr         1.0
            nChrp    "==1"
            nChrn    "==1"
            nNet     "==0"
          }
          .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=4"   # two radiative photons plus the two photons from pi0 -> gamma gamma
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :pion, against: [:kaon, :proton]   # only the lower-momentum track must be identified
            npip ">=1"
            npim ">=1"
          }
          # Reconstruct the pi0 from a photon pair: the two photons from the pi0 decay are
          # chosen as those minimising sqrt( ((M_gg - M_pi0)/sigma_pi0)^2 +
          # ((M_pi+pi-gammagamma - M_omega)/sigma_omega)^2 ), with sigma_pi0 ~ 7 MeV/c^2 and
          # sigma_omega ~ 6 MeV/c^2. This selection is realised by the chi2-minimising
          # combination choice of the Kalman pi0 reconstruction.
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
          }
          # 5C kinematic fit to psi' -> gamma_l gamma_h omega with the pi0 mass constrained.
          .kinematic_fit([:gamma, :gamma, :pi0, :pip, :pim]) {
            nominal
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            constrain_four_momentum
            chi2_cut 200   # loose BOSS cut; the paper requires chi2_5C <= 100 (applied in ROOT)
          }

alg_omega.note(:omega_signal_and_sideband,
               "omega signal region |M(pi+pi-pi0) - M_omega| <= 0.035 GeV/c^2; sidebands " \
               "0.68 <= M(pi+pi-pi0) <= 0.71 GeV/c^2 and 0.85 <= M(pi+pi-pi0) <= 0.88 GeV/c^2. " \
               "Windows applied in ROOT.")
         .note(:background_veto,
               "Background from multi-photon hadronic psi' decays suppressed by requiring " \
               "|M(gamma_l gamma_h) - M_eta| >= 25 MeV/c^2 and " \
               "|M(gamma_l gamma_h) - M_pi0| >= 15 MeV/c^2. The psi' -> gamma_h eta', " \
               "eta' -> gamma_l V background is suppressed by |M(gamma_l V) - M_eta'| > 15 MeV/c^2. " \
               "Applied in ROOT after the 5C fit.")
         .note(:mass_spectrum_fit,
               "The gamma_h omega invariant mass distribution is fitted with MC signal shapes for " \
               "the three chi_cJ resonances plus a background composed of the omega mass sideband " \
               "distribution and a 2nd-order polynomial. Fitted yields: 5 +- 11 (chi_c0, " \
               "<12.9e-6), 136 +- 14 (chi_c1, >10 sigma), 1 +- 6 (chi_c2, <6.1e-6). " \
               "Efficiencies 18.6% / 22.7% / 19.2%.")
         .note(:polarization_fit,
               "Transverse polarization fraction f_T from a likelihood fit to the cos(Theta) " \
               "distribution, Theta being the angle between the omega flight direction in the " \
               "chi_c1 rest frame and the normal to the omega decay plane in the omega rest " \
               "frame; f_T = 0.247 +0.090+0.044 -0.087-0.026 for chi_c1 -> gamma omega. ROOT-level.")
         .note(:systematics,
               "Tracking 4.0%, PID 2.0%, photon detection 4.0%, 4C fit 0.7%, 5C fit 3.1%, " \
               "selection-efficiency uncertainty 1.4%, background shape, binning, fit range and " \
               "sideband regions, signal shape, N_psi' and PDG branching fractions. Total ~9.4%.")

alg_omega.with_decay_card(decay_card_gammaomega).apply(sel_omega)
root_files_omega = alg_omega.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_gammaomega])
