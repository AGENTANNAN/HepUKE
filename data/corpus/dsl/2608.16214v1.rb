### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # (10087 +/- 44) x 10^6 J/psi events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

psip_data  = DatasetManager.real_data.find("709_3686")     # (2712 +/- 14) x 10^6 psi(3686) events
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ------------------------------------------------------------------------------
# Decay cards
# ------------------------------------------------------------------------------
decay_card_jpsi = <<~DECAYCARD
    Decay J/psi
    1.0000  Sigma0  anti-Sigma0  eta                    PHSP;
    Enddecay

    Decay Sigma0
    1.0000  gamma   Lambda0                             PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000  gamma   anti-Lambda0                        PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+   pi-                                    HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-   pi+                               HypWK;
    Enddecay

    Decay eta
    1.0000  gamma  gamma                                PHSP;
    Enddecay

    End
DECAYCARD

decay_card_psip = <<~DECAYCARD
    Decay psi(2S)
    1.0000  Sigma0  anti-Sigma0  eta                    PHSP;
    Enddecay

    Decay Sigma0
    1.0000  gamma   Lambda0                             PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000  gamma   anti-Lambda0                        PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+   pi-                                    HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-   pi+                               HypWK;
    Enddecay

    Decay eta
    1.0000  gamma  gamma                                PHSP;
    Enddecay

    End
DECAYCARD

exMC_jpsi_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_Sigma0Sigmabar_eta"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_jpsi
  config.cross_section   = :default
end

exMC_psip_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_Sigma0Sigmabar_eta"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_psip
  config.cross_section   = :default
end

# ==============================================================================
# Common event selection template (same track/photon/PID/Lambda/kfit chain).
# Different chi2_4C cut per resonance is applied in ROOT.
# ==============================================================================
def build_selection
  Selection.new
    .select_track {
       cos_theta 0.93
       Vz        20.0
       Vr        10.0
       nChrp     ">=2"
       nChrn     ">=2"
     }
    .select_photon {
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       tdc_emc_start     0
       tdc_emc_end       14
       nGam              ">=4"     # >= 4 good photons for Sigma0(->gamma Lambda) x2 + eta(->gg)
     }
    .pid(method: :probability) {
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]  # p and pbar; highest L(p)
       identify :pion,   against: [:kaon, :proton]
       nprp ">=1"
       nprm ">=1"
       npip ">=1"
       npim ">=1"
     }
    # Secondary vertex fits to reconstruct Lambda and Lambda_bar
    .secondary_vertex_fit([:prp, :pim]) {
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
    .secondary_vertex_fit([:prm, :pip]) {
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
    # Nominal 4C kinematic fit: J/psi (or psi') -> Lambda Lambdabar gamma gamma gamma gamma
    .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma, :gamma]) {
       nominal
       constrain_four_momentum
       chi2_cut 200                # loose in BOSS; tight cut (60 for J/psi, 40 for psi') in ROOT
     }
    # Alternative hypothesis: one fewer photon (for chi^2_ggg background rejection)
    .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma]) {
       constrain_four_momentum
     }
end

# ==============================================================================
# Algorithm for J/psi -> Sigma0 Sigmabar0 eta
# ==============================================================================
alg_jpsi = Algorithm.new("Sigma0SigmabarEtaJpsi")
alg_jpsi.set_header(["Sigma0SigmabarEtaJpsiAlg/Sigma0SigmabarEtaJpsi.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_jpsi = build_selection

alg_jpsi
  .note(:pi0_veto,
        "Photon-pair invariant mass veto |M(gg) - m(pi0)| > 20 MeV/c^2 (~3 sigma) " \
        "applied in ROOT to suppress pi0-containing backgrounds.")
  .note(:sigma_eta_selection,
        "Sigma0, Sigmabar0, eta candidates picked by minimising " \
        "chi^2 = sum ((M_gL - M_Sigma)/sigma_Sigma)^2 + ((M_gg - M_eta)/sigma_eta)^2 " \
        "with |M(gL) - m(Sigma0)| < 15 MeV/c^2, |M(g Lambdabar) - m(Sigmabar0)| < 15 " \
        "MeV/c^2 and |M(gg) - m(eta)| < 25 MeV/c^2, applied in ROOT.")
  .note(:ggg_hypothesis_veto,
        "Alternative 4C fit under J/psi -> Lambda Lambdabar gamma gamma gamma " \
        "used to reject events with chi^2_ggg < chi^2_4C (one-fewer-photon " \
        "background). Selection applied in ROOT.")

alg_jpsi.with_decay_card(decay_card_jpsi).apply(sel_jpsi)
root_files_jpsi = alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_signal])

# ==============================================================================
# Algorithm for psi(3686) -> Sigma0 Sigmabar0 eta
# ==============================================================================
alg_psip = Algorithm.new("Sigma0SigmabarEtaPsip")
alg_psip.set_header(["Sigma0SigmabarEtaPsipAlg/Sigma0SigmabarEtaPsip.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_psip = build_selection

alg_psip
  .note(:pi0_veto,
        "Photon-pair invariant mass veto |M(gg) - m(pi0)| > 20 MeV/c^2 (~3 sigma) " \
        "applied in ROOT to suppress pi0-containing backgrounds.")
  .note(:sigma_eta_selection,
        "Sigma0, Sigmabar0, eta candidates picked by minimising the joint " \
        "chi^2 with |M(gL)-m(Sigma0)|<15 MeV/c^2, |M(g Lambdabar)-m(Sigmabar0)|<15 " \
        "MeV/c^2, |M(gg)-m(eta)|<25 MeV/c^2 (ROOT).")
  .note(:jpsi_veto,
        "|M(Sigma0 Sigmabar0) - M(J/psi)| > 30 MeV/c^2 applied in ROOT to reject " \
        "psi(3686) -> eta J/psi, J/psi -> Sigma0 Sigmabar0 background.")

alg_psip.with_decay_card(decay_card_psip).apply(sel_psip)
root_files_psip = alg_psip.execute_on([psip_data, psip_incMC, exMC_psip_signal])
