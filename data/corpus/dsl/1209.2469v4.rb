# =============================================================================
# BESIII paper 1209.2469v4 — Search for invisible decays of eta and eta'
#   in J/psi -> phi eta and J/psi -> phi eta'
# arXiv:1209.2469v4
# Sample: (225.3 +- 2.8) x 10^6 J/psi events collected with BESIII
#
# Four event-selection chains:
#   1) J/psi -> phi eta,  phi -> K+ K-, eta  -> invisible   (recoil mass against phi)
#   2) J/psi -> phi eta', phi -> K+ K-, eta' -> invisible   (recoil mass against phi)
#   3) J/psi -> phi eta,  phi -> K+ K-, eta  -> gamma gamma (4C fit J/psi -> K+ K- gamma gamma)
#   4) J/psi -> phi eta', phi -> K+ K-, eta' -> gamma gamma (4C fit J/psi -> K+ K- gamma gamma)
# The eta(eta') -> gamma gamma channels are the normalisation channels from which
# the ratio B(eta(eta') -> invisible) / B(eta(eta') -> gamma gamma) is obtained.
# =============================================================================

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ---- Decay cards -------------------------------------------------------------

# Signal MC: J/psi -> phi eta, phi -> K+ K-, eta left UNDECAYED so that it is
# stable and escapes the detector (invisible decay).
decay_card_eta_inv = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta                  PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                    VSS;
    Enddecay

    End
DECAYCARD

# Signal MC: J/psi -> phi eta', phi -> K+ K-, eta' left undecayed (invisible).
decay_card_etap_inv = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta'                 PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                    VSS;
    Enddecay

    End
DECAYCARD

# Normalisation MC: J/psi -> phi eta, phi -> K+ K-, eta -> gamma gamma.
decay_card_eta_gg = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta                  PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                    VSS;
    Enddecay

    Decay eta
    1.0000 gamma gamma              PHSP;
    Enddecay

    End
DECAYCARD

# Normalisation MC: J/psi -> phi eta', phi -> K+ K-, eta' -> gamma gamma.
decay_card_etap_gg = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta'                 PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                    VSS;
    Enddecay

    Decay eta'
    1.0000 gamma gamma              PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive background (dominant for eta -> invisible):
# J/psi -> gamma eta_c, eta_c -> K+ pi- K_L0 (charge conjugate included in the fit).
decay_card_bkg_etac = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c              PHSP;
    Enddecay

    Decay eta_c
    1.0000 K+ pi- K_L0              PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive background (dominant for eta' -> invisible):
# J/psi -> phi K_L0 K_L0 (non-f_0(980) component).
decay_card_bkg_phiKLKL = <<~DECAYCARD
    Decay J/psi
    1.0000 phi K_L0 K_L0            PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                    VSS;
    Enddecay

    End
DECAYCARD

# Exclusive background (dominant for eta' -> invisible):
# J/psi -> phi f_0(980), f_0(980) -> K_L0 K_L0 (Flatte line shape in the fit).
decay_card_bkg_phif0 = <<~DECAYCARD
    Decay J/psi
    1.0000 phi f_0(980)             PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                    VSS;
    Enddecay

    Decay f_0(980)
    1.0000 K_L0 K_L0                PHSP;
    Enddecay

    End
DECAYCARD

# ---- Exclusive MC samples ----------------------------------------------------
exMC_eta_inv = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_phi_eta_invisible"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_eta_inv
  c.cross_section   = :default
end

exMC_etap_inv = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_phi_etap_invisible"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_etap_inv
  c.cross_section   = :default
end

exMC_eta_gg = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_phi_eta_gg"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_eta_gg
  c.cross_section   = :default
end

exMC_etap_gg = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_phi_etap_gg"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_etap_gg
  c.cross_section   = :default
end

exMC_bkg_etac = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_gamma_etac_KpiKL"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_bkg_etac
  c.cross_section   = :default
end

exMC_bkg_phiKLKL = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_phi_KLKL"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_bkg_phiKLKL
  c.cross_section   = :default
end

exMC_bkg_phif0 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_phi_f0980_KLKL"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_bkg_phif0
  c.cross_section   = :default
end

# =============================================================================
# ALGORITHM 1: J/psi -> phi eta, phi -> K+ K-, eta -> invisible
# The phi is tagged through K+ K-; the eta is undetected and inferred from the
# mass recoiling against the phi candidate.
# =============================================================================
alg_eta_inv = Algorithm.new("JpsiPhiEtaInvisible")
alg_eta_inv.set_header(["JpsiPhiEtaInvisibleAlg/JpsiPhiEtaInvisible.h"])
           .set_constant({"ECMS" => [:double, 3.097]})

sel_eta_inv = Selection.new
sel_eta_inv.select_track {
            # only the two charged tracks of the phi -> K+ K- candidate, net charge zero
            cos_theta 0.93
            Vz        10.0
            Vr         1.0
            nChrp     "==1"
            nChrn     "==1"
            nNet      "==0"
          }
          .select_photon {
            # EMC shower quality; no photon is required for the invisible channel,
            # the shower multiplicity in a cone around the recoil direction is vetoed
            # (see the notes below)
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion, :proton]   # ProbPID(K) > ProbPID(pi)
            nkp "==1"
            nkm "==1"
          }
          # Partial reconstruction: only the phi -> K+ K- tag is reconstructed; the
          # eta is an undetected particle inferred from the recoil four-momentum.
          # M^2_recoil(phi) = (p_lab - p_KK)^2 is the observable of the analysis.
          # rec_id_list: 0 => J/psi (top mother), 1 => phi, 2 => eta, 3 => K+, 4 => K-
          .partial_miss([2]) do               # rec_id 2 = eta (invisible / undetected)
            require_recoil_mass 0.495, 0.601  # +-3 sigma signal region around m_eta (sigma = 17.8 MeV)
          end

alg_eta_inv
  .note(:invisible_signal_mc,
        "the signal MC for J/psi -> phi eta, eta -> invisible is generated with the eta left " \
        "undecayed, so that it is stable and escapes the detector without leaving any signal.")
  .note(:shower_cone_veto,
        "the number of EMC showers N_shower (which could originate from a K_L or a photon) is " \
        "required to be zero inside a cone of 1.0 rad around the recoil direction against the phi " \
        "candidate; this rejects most eta decays into visible final states and most backgrounds " \
        "from multibody J/psi -> phi + anything decays. Applied in ROOT on the EMC shower list.")
  .note(:shower_energy_transition_region,
        "EMC showers are required to have E > 25 MeV in the barrel (|cos theta| < 0.80), " \
        "E > 50 MeV in the end-cap (0.86 < |cos theta| < 0.92) and E > 100 MeV in the barrel/" \
        "end-cap transition region (0.80 < |cos theta| < 0.86); showers must be isolated from all " \
        "charged tracks by more than 10 degrees.")
  .note(:recoil_direction_fiducial,
        "the polar angle of the recoil three-momentum against the phi candidate is required to " \
        "satisfy |cos(theta_recoil)| < 0.7 so that eta decay products are inside the detector " \
        "fiducial volume. Applied in ROOT.")
  .note(:phi_mass_window,
        "the phi tag requires 1.01 < m(K+ K-) < 1.03 GeV/c^2. Applied in ROOT.")
  .note(:recoil_mass_variable,
        "the mass recoiling against the phi candidate is computed with the four-momentum of the " \
        "incident beams in the lab frame, (M^recoil_phi)^2 = (p_lab - p_KK)^2 with " \
        "p_KK = p_phi = p_K+ + p_K-. The eta signal region is defined as +-3 sigma around the " \
        "known eta mass with a detector resolution sigma = 17.8 MeV/c^2 determined from MC.")
  .note(:background_estimation,
        "backgrounds are studied with more than 20 exclusive MC modes and divided into two " \
        "classes. Class I: J/psi -> phi eta with the eta decaying into visible final states that " \
        "are not detected by the EMC; expected 0.18 +- 0.02 events in the signal region. " \
        "Class II: J/psi decays without an eta (or without both eta and phi); the dominant " \
        "contribution is J/psi -> gamma eta_c, eta_c -> K+- pi-+ K_L0, where the soft radiative " \
        "photon is undetected or outside the 1 rad cone and the fast pion is misidentified as a " \
        "kaon; expected 0.8 +- 0.2 events. The eta_c -> K+- pi-+ K_L0 Dalitz-plot structure " \
        "uncertainty (phase-space vs. K_0*(1430) dominated) is estimated with a data-driven " \
        "reweighting of the simulation.")
  .note(:signal_extraction,
        "after all selection criteria, one event survives in the eta signal region where " \
        "1.0 +- 0.2 background events are expected. A 90% C.L. upper limit N_UL = 3.34 on the " \
        "number of J/psi -> phi eta, phi -> K+ K-, eta -> invisible events is obtained with the " \
        "POLE++ program using the Feldman-Cousins frequentist approach. ROOT level.")
  .note(:systematics,
        "relative systematic uncertainties on the ratio B(eta -> invisible)/B(eta -> gamma gamma): " \
        "N_shower requirement 0.3%, phi mass window 1.5%, J/psi -> gamma eta_c, eta_c -> K K pi " \
        "background 1.2%, 4C fit for eta -> gamma gamma 0.8%, photon detection 2.0%, signal shape " \
        "for eta -> gamma gamma 0.1%, background shape for eta -> gamma gamma 0.1%; " \
        "total systematic error 2.8%.")

alg_eta_inv.with_decay_card(decay_card_eta_inv).apply(sel_eta_inv)
root_files_eta_inv = alg_eta_inv.execute_on([jpsi_data, jpsi_incMC, exMC_eta_inv, exMC_bkg_etac])

# =============================================================================
# ALGORITHM 2: J/psi -> phi eta', phi -> K+ K-, eta' -> invisible
# =============================================================================
alg_etap_inv = Algorithm.new("JpsiPhiEtapInvisible")
alg_etap_inv.set_header(["JpsiPhiEtapInvisibleAlg/JpsiPhiEtapInvisible.h"])
            .set_constant({"ECMS" => [:double, 3.097]})

sel_etap_inv = Selection.new
sel_etap_inv.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr         1.0
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
           }
           .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             nkp "==1"
             nkm "==1"
           }
           # rec_id_list: 0 => J/psi (top mother), 1 => phi, 2 => eta', 3 => K+, 4 => K-
           .partial_miss([2]) do               # rec_id 2 = eta' (invisible / undetected)
             require_recoil_mass 0.80, 1.20    # M^recoil(phi) range of the unbinned extended ML fit
           end

alg_etap_inv
  .note(:invisible_signal_mc,
        "the signal MC for J/psi -> phi eta', eta' -> invisible is generated with the eta' left " \
        "undecayed, so that it is stable and escapes the detector without leaving any signal.")
  .note(:shower_cone_veto,
        "N_shower = 0 is required inside a cone of 1.0 rad around the recoil direction against " \
        "the phi candidate, rejecting eta' decays into visible final states and multibody " \
        "J/psi -> phi + anything backgrounds. Applied in ROOT.")
  .note(:shower_energy_transition_region,
        "EMC showers require E > 25 MeV (barrel |cos theta| < 0.80), E > 50 MeV (end-cap " \
        "0.86 < |cos theta| < 0.92) and E > 100 MeV in the transition region; showers must be " \
        "more than 10 degrees away from all charged tracks.")
  .note(:recoil_direction_fiducial,
        "the recoil direction against the phi candidate is required to satisfy " \
        "|cos(theta_recoil)| < 0.7. Applied in ROOT.")
  .note(:phi_mass_window,
        "the phi tag requires 1.01 < m(K+ K-) < 1.03 GeV/c^2. Applied in ROOT.")
  .note(:recoil_mass_variable,
        "(M^recoil_phi)^2 = (p_lab - p_KK)^2 with p_KK = p_K+ + p_K-. The eta' signal region is " \
        "within +-3 sigma of the known eta' mass, with a detector resolution sigma = 9.3 MeV/c^2 " \
        "determined from MC.")
  .note(:background_estimation,
        "Class I: J/psi -> phi eta' with the eta' decaying into undetected visible final states; " \
        "expected 1.0 +- 0.2 events. Class II: the dominant background is J/psi -> phi K_L K_L and " \
        "J/psi -> phi f_0(980), f_0(980) -> K_L K_L; expected 9.4 +- 1.7 events. The f_0(980) line " \
        "shape is parameterised with a Flatte form whose parameters are taken from BESII " \
        "measurements of J/psi -> phi pi+ pi- and phi K+ K-.")
  .note(:signal_extraction,
        "an unbinned extended maximum-likelihood fit to the M^recoil_phi distribution in the range " \
        "0.8 < M^recoil_phi < 1.2 GeV/c^2 is performed. The signal shape is fixed to the smoothed " \
        "MC histogram of a nearly background-free J/psi -> phi eta', eta' -> pi+ pi- eta, " \
        "eta -> gamma gamma control sample (purity > 98.5%), the J/psi -> phi f_0(980), " \
        "f_0(980) -> K_L K_L background shape and yield are fixed/floated from MC, and the " \
        "remaining J/psi -> phi K_L K_L background is modelled with a Chebychev polynomial. The " \
        "fitted signal yield N_sig = 2.3 +- 4.3 is consistent with zero; the 90% C.L. upper limit " \
        "on the number of signal events, obtained by integrating the normalised likelihood over " \
        "positive values, is N_UL = 10.1. ROOT level.")
  .note(:systematics,
        "relative systematic uncertainties on the ratio B(eta' -> invisible)/B(eta' -> gamma gamma): " \
        "N_shower requirement 0.3%, phi mass window 1.5%, J/psi -> phi f_0(980) background shape " \
        "1.0%, J/psi -> phi K_L K_L background shape 2.9%, 4C fit for eta' -> gamma gamma 0.8%, " \
        "photon detection 2.0%, signal shape for eta' -> gamma gamma 1.0%, background shape for " \
        "eta' -> gamma gamma 0.6%; total systematic error 4.1%.")

alg_etap_inv.with_decay_card(decay_card_etap_inv).apply(sel_etap_inv)
root_files_etap_inv = alg_etap_inv.execute_on([jpsi_data, jpsi_incMC, exMC_etap_inv,
                                               exMC_bkg_phiKLKL, exMC_bkg_phif0])

# =============================================================================
# ALGORITHM 3: J/psi -> phi eta, phi -> K+ K-, eta -> gamma gamma
# Normalisation channel; 4C kinematic fit under the J/psi -> K+ K- gamma gamma hypothesis.
# =============================================================================
alg_eta_gg = Algorithm.new("JpsiPhiEtaGG")
alg_eta_gg.set_header(["JpsiPhiEtaGGAlg/JpsiPhiEtaGG.h"])
          .set_constant({"ECMS" => [:double, 3.097]})

sel_eta_gg = Selection.new
sel_eta_gg.select_track {
            # same charged-track criteria as the invisible channel
            cos_theta 0.93
            Vz        10.0
            Vr         1.0
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
            nGam              ">=2"   # at least two good photons for eta -> gamma gamma
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion, :proton]
            nkp "==1"
            nkm "==1"
          }
          # 4C kinematic fit under the J/psi -> K+ K- gamma gamma hypothesis; if more than
          # two photons are present the fit is repeated over all permutations and the
          # combination with the best fit is retained (handled automatically).
          .kinematic_fit([:kp, :km, :gamma, :gamma]) {
            nominal
            constrain_four_momentum
            chi2_cut 200   # paper: chi2(K+K-gammagamma) < 90 for the eta case (applied in ROOT)
          }

alg_eta_gg
  .note(:kinematic_fit_chi2,
        "the J/psi -> K+ K- gamma gamma 4C kinematic fit chi2 is required to be less than 90 for " \
        "the eta case; the optimal cut is applied in ROOT.")
  .note(:mass_windows,
        "the unbinned extended ML fit is performed over 0.99 < m(K+ K-) < 1.10 GeV/c^2 and " \
        "0.35 < m(gamma gamma) < 0.75 GeV/c^2. ROOT level.")
  .note(:signal_and_background_shapes,
        "the phi signal shape is a relativistic Breit-Wigner convolved with a Gaussian detector " \
        "resolution (phi width fixed to the PDG value, mass floated); the eta signal shape is a " \
        "Crystal Ball function with floated parameters. Backgrounds are divided into three " \
        "categories: non-phi-eta-peaking (J/psi -> gamma pi0 K+ K-, one photon missing), " \
        "non-phi-peaking (J/psi -> K+ K- eta), and non-eta-peaking (J/psi -> phi gamma gamma and " \
        "phi pi0 pi0). The non-phi-peaking m(K+K-) shape is parameterised as " \
        "B(m_KK) = (m_KK - 2 m_K)^a exp(-b m_KK - c m_KK^2); the non-eta-peaking m(gamma gamma) " \
        "shape is a second-order Chebychev polynomial. The yields and the ratio " \
        "B(J/psi -> phi eta)/B(J/psi -> phi eta') are extracted from a simultaneous fit. " \
        "ROOT level.")
  .note(:background_shapes,
        "background shapes in m(K+ K-) and m(gamma gamma) are floated in the fit to data. " \
        "The background PDFs are the non-phi-eta-peaking, non-phi-peaking and non-eta-peaking " \
        "components listed above. ROOT level.")
  .note(:systematics,
        "relative systematic uncertainties on the eta -> gamma gamma normalisation: photon " \
        "detection 1% per photon (2.0% for the two-photon final state), 4C fit 0.4%, phi mass " \
        "resolution 0.1%, non-phi-eta-peaking background shape 0.1%; the remaining uncertainties " \
        "on the ratios (N_shower requirement, phi mass window, background shapes) are quoted in " \
        "Table II of the paper. ROOT level.")

alg_eta_gg.with_decay_card(decay_card_eta_gg).apply(sel_eta_gg)
root_files_eta_gg = alg_eta_gg.execute_on([jpsi_data, jpsi_incMC, exMC_eta_gg])

# =============================================================================
# ALGORITHM 4: J/psi -> phi eta', phi -> K+ K-, eta' -> gamma gamma
# Normalisation channel; 4C kinematic fit under the J/psi -> K+ K- gamma gamma hypothesis.
# =============================================================================
alg_etap_gg = Algorithm.new("JpsiPhiEtapGG")
alg_etap_gg.set_header(["JpsiPhiEtapGGAlg/JpsiPhiEtapGG.h"])
           .set_constant({"ECMS" => [:double, 3.097]})

sel_etap_gg = Selection.new
sel_etap_gg.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr         1.0
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
             nGam              ">=2"   # at least two good photons for eta' -> gamma gamma
           }
           .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             nkp "==1"
             nkm "==1"
           }
           .kinematic_fit([:kp, :km, :gamma, :gamma]) {
             nominal
             constrain_four_momentum
             chi2_cut 200   # paper: chi2(K+K-gammagamma) < 40 for the eta' case (applied in ROOT)
           }

alg_etap_gg
  .note(:kinematic_fit_chi2,
        "the J/psi -> K+ K- gamma gamma 4C kinematic fit chi2 is required to be less than 40 for " \
        "the eta' case; the optimal cut is applied in ROOT.")
  .note(:mass_windows,
        "the unbinned extended ML fit is performed over 0.99 < m(K+ K-) < 1.10 GeV/c^2 and " \
        "0.75 < m(gamma gamma) < 1.15 GeV/c^2. ROOT level.")
  .note(:signal_and_background_shapes,
        "the phi signal shape is a relativistic Breit-Wigner convolved with a Gaussian resolution " \
        "(width fixed to PDG, mass floated); the eta' signal shape is a Crystal Ball function with " \
        "floated parameters. The backgrounds are the non-phi-eta'-peaking (J/psi -> gamma pi0 K+ K-), " \
        "non-phi-peaking (J/psi -> K+ K- eta') and non-eta'-peaking (J/psi -> phi gamma gamma and " \
        "phi pi0 pi0) components; the non-phi-peaking m(K+ K-) shape uses " \
        "B(m_KK) = (m_KK - 2 m_K)^a exp(-b m_KK - c m_KK^2) and the non-eta'-peaking m(gamma gamma) " \
        "shape a second-order Chebychev polynomial. Yields are extracted from the simultaneous fit. " \
        "ROOT level.")
  .note(:systematics,
        "relative systematic uncertainties on the eta' -> gamma gamma normalisation: photon " \
        "detection 1% per photon (2.0%), 4C fit 0.8%, phi mass resolution 1.0%, " \
        "non-phi-eta'-peaking background shape 0.6%; the total systematic error on the ratio is " \
        "4.1%. ROOT level.")

alg_etap_gg.with_decay_card(decay_card_etap_gg).apply(sel_etap_gg)
root_files_etap_gg = alg_etap_gg.execute_on([jpsi_data, jpsi_incMC, exMC_etap_gg])
