### Dataset description ###
# ψ(3686) real data and its corresponding inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# --- Decay cards (EvtGen) ---
# The four radiative signals all share the final state γ p p̄ K⁺K⁻.
# The radiative transition ψ(3686) → γ X is generated with VSP_PWAVE;
# the subsequent X → p p̄ K⁺K⁻ is generated with PHSP.
decay_card_etac2S = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c(2S) VSP_PWAVE;
    Enddecay

    Decay eta_c(2S)
    1.000 p+ anti-p- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 VSP_PWAVE;
    Enddecay

    Decay chi_c0
    1.000 p+ anti-p- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 VSP_PWAVE;
    Enddecay

    Decay chi_c1
    1.000 p+ anti-p- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 VSP_PWAVE;
    Enddecay

    Decay chi_c2
    1.000 p+ anti-p- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Non-radiative background: ψ(3686) → p p̄ K⁺K⁻
decay_card_bkg = <<~DECAYCARD
    Decay psi(2S)
    1.000 p+ anti-p- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC samples ---
# 500k events for each of the four radiative signal modes
exMC_etac2S = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gam_etac2S_ppKK"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_etac2S
  config.cross_section   = :default
end

exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gam_chic0_ppKK"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic0
  config.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gam_chic1_ppKK"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic1
  config.cross_section   = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gam_chic2_ppKK"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic2
  config.cross_section   = :default
end

# 200k events for the non-radiative ψ(3686) → p p̄ K⁺K⁻ background
exMC_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_ppKK_bkg"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_bkg
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# All four radiative signals share the final state γ p p̄ K⁺K⁻ and the same
# selection criteria, so a single Algorithm object covers all of them.
alg_name = "PsipGammaPPKK"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93     # |cosθ| < 0.93
                  Vz        10.0     # |Vz| < 10 cm
                  Vr        1.0      # Vr < 1 cm in the transverse plane
                  nChrp    "==2"     # exactly two positive tracks
                  nChrn    "==2"     # exactly two negative tracks
                  nNet     "==0"     # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0      # TDC window start (=> 0 ns)
                  tdc_emc_end       14     # TDC window end   (=> 700 ns)
                  energyThreshold_b 0.025  # E > 25 MeV in the barrel
                  energyThreshold_e 0.050  # E > 50 MeV in the endcap
                  nGam              ">=1"  # at least one good photon
                }
               # Step 1: identify one proton and one anti-proton
               # (against kaon and pion hypotheses)
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # p+ and p̄ at once
                  nprp ">=1"
                  nprm ">=1"
                }
               # ... then remove the proton / anti-proton tracks
               .remove([:prp <= :chrgp, :prm <= :chrgn])
               # Step 2: identify one K⁺ and one K⁻ (against pion) from the remaining tracks
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :kaon, against: [:pion]  # K+ and K- at once
                  nkp ">=1"
                  nkm ">=1"
                }
               # 3C kinematic fit on γ p p̄ K⁺K⁻: a common vertex is fitted on the
               # four charged tracks and only the three-momentum is constrained,
               # leaving the photon energy floating. The fit is flagged nominal.
               .kinematic_fit([:gamma, :prp, :prm, :kp, :km]) {
                  nominal
                  vertex_fit([1, 2, 3, 4])        # common vertex for prp, prm, kp, km
                  constrain_three_momentum        # 3C constraint (photon energy floating)
                  chi2_cut 200
                }

# --- Inexpressible BOSS-side procedures preserved as notes ---
my_algorithm
  .note(:helix_correction, "helix-parameter track correction applied to all charged
    tracks before the 3C kinematic fit")
  .note(:fsr_correction, "final-state-radiation correction applied to the signal selection")
  .note(:fake_photon_resampling, "fake-photon resampling used to model the background
    photon contribution in the γ p p̄ K⁺K⁻ final state")
  .note(:efficiency_curve, "GPR-based efficiency curve used to correct the reconstruction
    efficiency as a function of the kinematic variable")
  .note(:photon_angular_distribution, "signal MC generated with a 1 + alpha*cos^2(theta_gamma)
    photon angular distribution")
  .note(:lineshape_damping, "damping function applied in the signal lineshape")
  .note(:resolution_smearing, "MC-data resolution smearing applied to match the data resolution")

# The header/kinematic variables generated from this card (γ, p, p̄, K⁺, K⁻) are
# identical for all four radiative decay cards, since the final state is shared.
my_algorithm.with_decay_card(decay_card_etac2S).apply(event_selection)

root_files = my_algorithm.execute_on([psip_data, psip_incMC,
                                      exMC_etac2S, exMC_chic0, exMC_chic1, exMC_chic2,
                                      exMC_bkg])