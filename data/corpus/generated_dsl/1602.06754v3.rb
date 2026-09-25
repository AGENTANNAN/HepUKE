# =============================================================================
# BESIII: single-baryon-tag reactions  psi -> Xi- Xi+  and
#         psi -> Sigma(1385)- Sigma(1385)+  (charge conjugates included),
# on J/psi (3.097 GeV) and psi(3686) (3.686 GeV).
# The tag baryon is reconstructed as Xi-/Sigma(1385)- -> pi- Lambda, Lambda -> p pi-,
# and the anti-baryon is inferred from the pi-Lambda recoil mass (partial
# reconstruction; no kinematic fit is performed).
# =============================================================================

### Datasets ###
jpsi_data  = DatasetManager.real_data.find("708_3097")        # 223.7M J/psi events @ 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")     # inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")        # 106.4M psi(3686) events @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")     # inclusive MC

### Decay cards ###
# The anti-baryon is left as a stable final state: it is never reconstructed,
# it is inferred from the recoil mass of the pi-Lambda (tag-baryon) system.
decay_card_jpsi_xi = <<~DECAYCARD
    Decay J/psi
    1.0 Xi- anti-Xi+ PHSP;
    Enddecay
    Decay Xi-
    1.0 pi- Lambda0 HypWK;
    Enddecay
    Decay Lambda0
    1.0 p+ pi- HypWK;
    Enddecay
    End
DECAYCARD

decay_card_psip_xi = <<~DECAYCARD
    Decay psi(2S)
    1.0 Xi- anti-Xi+ PHSP;
    Enddecay
    Decay Xi-
    1.0 pi- Lambda0 HypWK;
    Enddecay
    Decay Lambda0
    1.0 p+ pi- HypWK;
    Enddecay
    End
DECAYCARD

decay_card_jpsi_sigma = <<~DECAYCARD
    Decay J/psi
    1.0 Sigma(1385)- anti-Sigma(1385)+ PHSP;
    Enddecay
    Decay Sigma(1385)-
    1.0 pi- Lambda0 PHSP;
    Enddecay
    Decay Lambda0
    1.0 p+ pi- HypWK;
    Enddecay
    End
DECAYCARD

decay_card_psip_sigma = <<~DECAYCARD
    Decay psi(2S)
    1.0 Sigma(1385)- anti-Sigma(1385)+ PHSP;
    Enddecay
    Decay Sigma(1385)-
    1.0 pi- Lambda0 PHSP;
    Enddecay
    Decay Lambda0
    1.0 p+ pi- HypWK;
    Enddecay
    End
DECAYCARD

### Exclusive MC: 1M events for each of the four modes ###
exMC_jpsi_xi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3097_xi_xi"
    config.related_dataset = jpsi_data
    config.events          = 1_000_000
    config.decay_card      = decay_card_jpsi_xi
    config.cross_section   = :default
end

exMC_psip_xi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3686_xi_xi"
    config.related_dataset = psip_data
    config.events          = 1_000_000
    config.decay_card      = decay_card_psip_xi
    config.cross_section   = :default
end

exMC_jpsi_sigma = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3097_sigma_sigma"
    config.related_dataset = jpsi_data
    config.events          = 1_000_000
    config.decay_card      = decay_card_jpsi_sigma
    config.cross_section   = :default
end

exMC_psip_sigma = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3686_sigma_sigma"
    config.related_dataset = psip_data
    config.events          = 1_000_000
    config.decay_card      = decay_card_psip_sigma
    config.cross_section   = :default
end

# =============================================================================
# Mode 1: J/psi -> Xi- Xi+   (tag Xi- -> pi- Lambda, Lambda -> p pi-)
# =============================================================================
alg_jpsi_xi = Algorithm.new("JpsiXiRecoil")
alg_jpsi_xi.set_header(["JpsiXiRecoilAlg/JpsiXiRecoil.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .note(:vertex_fit_chi2,
                 "Lambda (p pi-) secondary vertex fit required chi2 < 500")
           .note(:angular_distribution,
                 "exclusive MC generation must incorporate the measured alpha angular "
                 "distribution of J/psi -> Xi- Xi+")

sel_jpsi_xi = Selection.new
    .select_track {                               # charged track quality + multiplicity
        cos_theta 0.93                            # |cos(theta)| < 0.93
        Vz        10.0                            # |Vz| < 10 cm
        Vr        1.0                             # Vr < 1 cm
        nChrp     ">=1"                           # >=1 positive track (p) ...
        nChrn     ">=1"                           # ... >=1 negative track (pi-)
        nTot      ">=3"                           # union of (1p,2pi-) and its charge conjugate (2pi+,1pbar)
    }
    .pid(method: :probability) {                  # PID, probability method
        prob_cut 0.001                            # prob > 0.001
        identify :pion,   against: [:kaon, :proton]   # pi vs K, p  (tag side has 2 pi-)
        identify :proton, against: [:kaon, :pion]     # p  vs K, pi
        npim ">=2"                                # >=2 pi
        nprp ">=1"                                # >=1 p
    }
    .secondary_vertex_fit([:prp, :pim]) {         # p pi- common vertex -> Lambda
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        invariant_mass_of(:prp, :pim).within(1.1097, 1.1217)  # |M(p pi-) - M_Lambda| < 6 MeV/c^2
        remove_used_particle_from_candidate_list  # do not reuse the Lambda daughters
    }
    # No kinematic fit: a partial reconstruction tags pi- Lambda and misses anti-Xi+
    # (decay-card recID 2), whose mass is read from the pi-Lambda recoil mass.
    .partial_miss([2]) {
        best_combination_by_mass :Xi-, 1.3217     # pick the pi-Lambda combo closest to M_Xi
        require_recoil_mass 1.312, 1.332          # pi-Lambda recoil mass window (Xi)
    }

alg_jpsi_xi.with_decay_card(decay_card_jpsi_xi).apply(sel_jpsi_xi)
alg_jpsi_xi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_xi])

# =============================================================================
# Mode 2: psi(3686) -> Xi- Xi+   (tag Xi- -> pi- Lambda, Lambda -> p pi-)
# =============================================================================
alg_psip_xi = Algorithm.new("PsipXiRecoil")
alg_psip_xi.set_header(["PsipXiRecoilAlg/PsipXiRecoil.h"])
           .set_constant({"ECMS" => [:double, 3.686]})
           .note(:vertex_fit_chi2,
                 "Lambda (p pi-) secondary vertex fit required chi2 < 500")
           .note(:angular_distribution,
                 "exclusive MC generation must incorporate the measured alpha angular "
                 "distribution of psi(3686) -> Xi- Xi+")
           .note(:background_veto,
                 "offline veto |M_recoil(pi+ pi-) - M_J/psi| > 0.005 GeV/c^2 to suppress "
                 "psi(3686) -> pi+ pi- J/psi; applied at the ROOT analysis stage")

sel_psip_xi = Selection.new
    .select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=1"
        nChrn     ">=1"
        nTot      ">=3"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :pion,   against: [:kaon, :proton]
        identify :proton, against: [:kaon, :pion]
        npim ">=2"
        nprp ">=1"
    }
    .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        invariant_mass_of(:prp, :pim).within(1.1097, 1.1217)
        remove_used_particle_from_candidate_list
    }
    .partial_miss([2]) {
        best_combination_by_mass :Xi-, 1.3217
        require_recoil_mass 1.308, 1.338          # pi-Lambda recoil mass window (psi(3686) Xi)
    }

alg_psip_xi.with_decay_card(decay_card_psip_xi).apply(sel_psip_xi)
alg_psip_xi.execute_on([psip_data, psip_incMC, exMC_psip_xi])

# =============================================================================
# Mode 3: J/psi -> Sigma(1385)- Sigma(1385)+
#         (tag Sigma(1385)- -> pi- Lambda, Lambda -> p pi-)
# =============================================================================
alg_jpsi_sigma = Algorithm.new("JpsiSigmaRecoil")
alg_jpsi_sigma.set_header(["JpsiSigmaRecoilAlg/JpsiSigmaRecoil.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .note(:vertex_fit_chi2,
                    "Lambda (p pi-) secondary vertex fit required chi2 < 500")
              .note(:angular_distribution,
                    "exclusive MC generation must incorporate the measured alpha angular "
                    "distribution of J/psi -> Sigma(1385)- Sigma(1385)+")

sel_jpsi_sigma = Selection.new
    .select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=1"
        nChrn     ">=1"
        nTot      ">=3"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :pion,   against: [:kaon, :proton]
        identify :proton, against: [:kaon, :pion]
        npim ">=2"
        nprp ">=1"
    }
    .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        invariant_mass_of(:prp, :pim).within(1.1097, 1.1217)
        remove_used_particle_from_candidate_list
    }
    .partial_miss([2]) {
        best_combination_by_mass :"Sigma(1385)-", 1.3872   # closest to M_Sigma(1385)
        require_recoil_mass 1.350, 1.420                   # pi-Lambda recoil mass window (Sigma(1385))
    }

alg_jpsi_sigma.with_decay_card(decay_card_jpsi_sigma).apply(sel_jpsi_sigma)
alg_jpsi_sigma.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_sigma])

# =============================================================================
# Mode 4: psi(3686) -> Sigma(1385)- Sigma(1385)+
#         (tag Sigma(1385)- -> pi- Lambda, Lambda -> p pi-)
# =============================================================================
alg_psip_sigma = Algorithm.new("PsipSigmaRecoil")
alg_psip_sigma.set_header(["PsipSigmaRecoilAlg/PsipSigmaRecoil.h"])
              .set_constant({"ECMS" => [:double, 3.686]})
              .note(:vertex_fit_chi2,
                    "Lambda (p pi-) secondary vertex fit required chi2 < 500")
              .note(:angular_distribution,
                    "exclusive MC generation must incorporate the measured alpha angular "
                    "distribution of psi(3686) -> Sigma(1385)- Sigma(1385)+")
              .note(:background_veto,
                    "offline veto |M_recoil(pi+ pi-) - M_J/psi| > 0.005 GeV/c^2 to suppress "
                    "psi(3686) -> pi+ pi- J/psi; applied at the ROOT analysis stage")

sel_psip_sigma = Selection.new
    .select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=1"
        nChrn     ">=1"
        nTot      ">=3"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :pion,   against: [:kaon, :proton]
        identify :proton, against: [:kaon, :pion]
        npim ">=2"
        nprp ">=1"
    }
    .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        invariant_mass_of(:prp, :pim).within(1.1097, 1.1217)
        remove_used_particle_from_candidate_list
    }
    .partial_miss([2]) {
        best_combination_by_mass :"Sigma(1385)-", 1.3872
        require_recoil_mass 1.350, 1.420
    }

alg_psip_sigma.with_decay_card(decay_card_psip_sigma).apply(sel_psip_sigma)
alg_psip_sigma.execute_on([psip_data, psip_incMC, exMC_psip_sigma])