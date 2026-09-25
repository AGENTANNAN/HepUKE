# BESIII: Λc+ → p η' (SCS decay) using single-tag method with DNN
# Also measures normalization channel Λc+ → p ω (ω → π+π-π0)
# η' reconstructed via η' → π+π-γ
# Uses 4.5 fb^-1 collected at √s ∈ [4.600, 4.699] GeV.

### Datasets — Λc+ Λc- threshold scan (2016-2019) ###
data_4600  = DatasetManager.real_data.find("703_4600")
data_4612  = DatasetManager.real_data.find("706_4610")
data_4626  = DatasetManager.real_data.find("706_4620")
data_4640  = DatasetManager.real_data.find("706_4640")
data_4660  = DatasetManager.real_data.find("706_4660")
data_4680  = DatasetManager.real_data.find("706_4680")
data_4700  = DatasetManager.real_data.find("706_4700")

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4612 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4626 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

data_all   = [data_4600, data_4612, data_4626, data_4640, data_4660, data_4680, data_4700]
incMC_all  = [incMC_4600, incMC_4612, incMC_4626, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

### Decay cards ###
# Signal: Λc+ → p η', η' → π+π-γ
decay_card_petap = <<~DECAYCARD
  Decay Lambda_c+
  1.0000  p+  eta'                             PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000  anti-p-  eta'                        PHSP;
  Enddecay

  Decay eta'
  1.0000  pi+  pi-  gamma                      ETA_DALITZ;
  Enddecay

  End
DECAYCARD

# Normalization: Λc+ → p ω, ω → π+π-π0
decay_card_pomega = <<~DECAYCARD
  Decay Lambda_c+
  1.0000  p+  omega                            PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000  anti-p-  omega                       PHSP;
  Enddecay

  Decay omega
  1.0000  pi+  pi-  pi0                        OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.0000  gamma  gamma                         PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC ###
exMC_petap = DatasetManager.create_exclusive_mc_for(data_all) do |c|
  c.sample_name     = "Lcp_p_etap"
  c.events          = 200000
  c.decay_card      = decay_card_petap
  c.cross_section   = :default
end
exMC_pomega = DatasetManager.create_exclusive_mc_for(data_all) do |c|
  c.sample_name     = "Lcp_p_omega"
  c.events          = 200000
  c.decay_card      = decay_card_pomega
  c.cross_section   = :default
end

######################################################################
# Algorithm 1: Signal — Λc+ → p η', η' → π+π-γ  (ST reconstruction)
######################################################################
alg_sig = Algorithm.new("LcpToPetap")
alg_sig.set_header(["LcpToPetapAlg/LcpToPetap.h"])

sel_sig = Selection.new
sel_sig.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp    ">=2"      # p and pi+
         nChrn    ">=1"      # pi-
         nTot     ">=3"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :proton, against: [:kaon, :pion]
         identify :pion,   against: [:kaon]
       }
       # η' reconstruction: M(π+π-γ) ∈ [0.90, 1.00] GeV/c^2  (mass-constrained)
       .kalman_kinematic_fit([:pip, :pim, :gamma]) {
         invariant_mass_of(:pip, :pim, :gamma).constrain_to_nominal_mass_of(:etap)
         chi2_cut 200
         netap ">=1"
       }
       # Λc+ mass-constraint on p η' (single-tag hypothesis).
       # This serves as the terminating kinematic-fit step.
       .kinematic_fit([:prp, :etap]) {
         nominal
         invariant_mass_of(:prp, :etap).constrain_to_nominal_mass_of(:"Lambda_c+")
         chi2_cut 200
       }

alg_sig.note(:eta_prime_mass_window,
             "M(π+π-γ) ∈ [0.90, 1.00] GeV/c^2 to reconstruct η'.")
       .note(:background_veto,
             "M(pπ-) required to be outside [1.10, 1.15] GeV/c^2 (Λ veto); " \
             "M(π+π-) required outside [0.48, 0.51] GeV/c^2 (K_S0 veto).")
       .note(:mbc_deltaE,
             "Λc+ candidate identified via M_BC = sqrt(E_beam^2 - |p_Λc|^2) " \
             "and ΔE = E_Λc - E_beam; retain combination with smallest |ΔE| " \
             "and require -65 MeV < ΔE < 60 MeV.")
       .note(:dnn_classifier,
             "Deep-learning classifier based on Particle Transformer (ParT) " \
             "with 10-model ensemble; three-class output (signal, ΛcΛc-bar " \
             "non-signal, hadronic). Composite discriminator " \
             "S = score_sig · (1 - score_hadronic) > 0.98 applied for final " \
             "event selection. Trained on ~1.5M events per channel; " \
             "background reduced by >2 orders of magnitude retaining ~40% signal.")

alg_sig.with_decay_card(decay_card_petap).apply(sel_sig)
alg_sig.execute_on(data_all + incMC_all + exMC_petap)

######################################################################
# Algorithm 2: Normalization — Λc+ → p ω, ω → π+π-π0
######################################################################
alg_norm = Algorithm.new("LcpToPomega")
alg_norm.set_header(["LcpToPomegaAlg/LcpToPomega.h"])

sel_norm = Selection.new
sel_norm.select_track {
          cos_theta 0.93
          Vz        10.0
          Vr        1.0
          nChrp    ">=2"       # p and pi+
          nChrn    ">=1"       # pi-
          nTot     ">=3"
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
          identify :proton, against: [:kaon, :pion]
          identify :pion,   against: [:kaon]
        }
        # π0 reconstruction: γγ mass window (0.115, 0.150) GeV/c^2, 1C mass fit
        .kalman_kinematic_fit([:gamma, :gamma]) {
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 200
          npi0 ">=1"
        }
        # Λc+ mass-constraint on p (π+π-π0) — nominal terminating kinematic fit
        .kinematic_fit([:prp, :pip, :pim, :pi0]) {
          nominal
          invariant_mass_of(:pip, :pim, :pi0).constrain_to_nominal_mass_of(:omega)
          invariant_mass_of(:prp, :pip, :pim, :pi0).constrain_to_nominal_mass_of(:"Lambda_c+")
          chi2_cut 200
        }

alg_norm.note(:omega_mass_window,
              "M(π+π-π0) ∈ [0.73, 0.83] GeV/c^2 to reconstruct ω.")
        .note(:background_veto,
              "M(pπ-) required outside [1.10, 1.15] GeV/c^2 (Λ veto); " \
              "M(π+π-) required outside [0.48, 0.51] GeV/c^2 (K_S0 veto); " \
              "M(pπ0) required outside [1.17, 1.20] GeV/c^2 (Σ+ veto).")
        .note(:mbc_deltaE,
              "Λc+ candidate identified via M_BC and ΔE; retain smallest " \
              "|ΔE| combination and require -40 MeV < ΔE < 30 MeV.")
        .note(:dnn_classifier,
              "Same DNN (Particle Transformer + ensemble of 10) as in signal " \
              "channel, jointly trained on both modes for consistency; " \
              "S = score_sig · (1 - score_hadronic) > 0.98.")

alg_norm.with_decay_card(decay_card_pomega).apply(sel_norm)
alg_norm.execute_on(data_all + incMC_all + exMC_pomega)
