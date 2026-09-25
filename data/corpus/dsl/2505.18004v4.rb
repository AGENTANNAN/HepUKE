# =====================================================================
# BESIII: Measurement of BF of Λc+ → Σ+ η and Σ+ η'
#   [arXiv:2505.18004v4]
#
# Analysis of e+e- data at 7 energy points (4.600–4.699 GeV, 4.5 fb⁻¹).
# Single-tag method: fully reconstruct Λc+ candidates from their
# decay products. Signal: Λc+ → Σ+ η (η→γγ, π+π-π0) and
# Λc+ → Σ+ η' (η'→π+π-η, η→γγ).
# Reference: Λc+ → Σ+ π0 and Λc+ → Σ+ ω (ω→π+π-π0).
# Ordinary analysis (Algorithm + Selection). Note: multiple independent
# decay channels each get their own Algorithm per Rule T1.
# =====================================================================

# Energy points (BOSS versions from BES3_dataset.md)
# 703_4600 = 4599.53 MeV
# 706_4610 = 4611.86 MeV
# 706_4620 = 4628.00 MeV
# 706_4640 = 4640.91 MeV
# 706_4660 = 4661.24 MeV
# 706_4680 = 4681.92 MeV
# 706_4700 = 4698.82 MeV

data_4600  = DatasetManager.real_data.find("703_4600")
data_4612  = DatasetManager.real_data.find("706_4610")
data_4628  = DatasetManager.real_data.find("706_4620")
data_4641  = DatasetManager.real_data.find("706_4640")
data_4661  = DatasetManager.real_data.find("706_4660")
data_4682  = DatasetManager.real_data.find("706_4680")
data_4699  = DatasetManager.real_data.find("706_4700")

all_data = [data_4600, data_4612, data_4628, data_4641, data_4661, data_4682, data_4699]

inc_mc_4600  = DatasetManager.inclusive_mc.find("703_4600")
inc_mc_4612  = DatasetManager.inclusive_mc.find("706_4610")
inc_mc_4628  = DatasetManager.inclusive_mc.find("706_4620")
inc_mc_4641  = DatasetManager.inclusive_mc.find("706_4640")
inc_mc_4661  = DatasetManager.inclusive_mc.find("706_4660")
inc_mc_4682  = DatasetManager.inclusive_mc.find("706_4680")
inc_mc_4699  = DatasetManager.inclusive_mc.find("706_4700")

all_inc_mc = [inc_mc_4600, inc_mc_4612, inc_mc_4628, inc_mc_4641, inc_mc_4661, inc_mc_4682, inc_mc_4699]

# Common sub-decay cards
sub_pi0_gg = <<~DECAYCARD
  Decay pi0
  1.0  gamma gamma  PHSP;
  Enddecay
DECAYCARD

sub_eta_gg = <<~DECAYCARD
  Decay eta
  1.0  gamma gamma  PHSP;
  Enddecay
DECAYCARD

sub_eta_3pi = <<~DECAYCARD
  Decay eta
  1.0  pi+ pi- pi0  PHSP;
  Enddecay
DECAYCARD

sub_omega_3pi = <<~DECAYCARD
  Decay omega
  1.0  pi+ pi- pi0  PHSP;
  Enddecay
DECAYCARD

sub_Sigma_p_pi0 = <<~DECAYCARD
  Decay Sigma+
  1.0  p pi0  PHSP;
  Enddecay
DECAYCARD

# =====================================================================
# Algorithm 1: Λc+ → Σ+ η (η → γγ), reference: Λc+ → Σ+ π0
# =====================================================================

decay_card_sig_eta_gg = <<~DECAYCARD
  NoDecay e+ e-

  Decay Lambda_c+
  1.0  Sigma+ eta  PHSP;
  Enddecay

  #{sub_Sigma_p_pi0}
  #{sub_eta_gg}
  #{sub_pi0_gg}
  End
DECAYCARD

exMC_sig_eta_gg = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name     = "Lc_to_Sigma_eta_gg"
  config.decay_card      = decay_card_sig_eta_gg
  config.events          = 200_000
  config.cross_section   = :default
end

alg1 = Algorithm.new("Lc2SigmaEtaGG")
alg1.set_header(["Lc2SigmaEtaGGAlg/Lc2SigmaEtaGG.h"])

event_selection_1 = Selection.new
event_selection_1.select_track {
                   cos_theta 0.93
                   Vz        10.0
                   Vr        1.0
                   nChrp     ">=1"
                   nChrn     ">=1"
                 }
                .select_photon {
                   tdc_emc_start     0
                   tdc_emc_end       14
                   angle_to_track    10.0
                   energyThreshold_b 0.025
                   energyThreshold_e 0.050
                   nGam              ">=4"    # 2 from π0 + 2 from η
                 }
                .pid(method: :probability) {
                   prob_cut 0.0
                   identify :prp, against: [:kp, :pip]   # proton for Σ+
                   identify :kp,  against: :pip
                   identify :pip,  against: :kp
                 }
                .kalman_kinematic_fit([:gamma, :gamma]) {
                   invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                   chi2_cut 100
                 }
                .kalman_kinematic_fit([:gamma, :gamma]) {
                   invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                   chi2_cut 100
                 }
                .kinematic_fit([:prp, :pi0, :eta]) {
                   nominal
                   constrain_four_momentum
                   invariant_mass_of(:prp, :pi0).constrain_to_nominal_mass_of(:"Sigma+")
                   chi2_cut 17     # optimized by FOM for Σ+η
                 }

alg1.with_decay_card(decay_card_sig_eta_gg)
    .note(:signal_region,
          "Signal region: M(pπ0) ∈ (1.174, 1.200) GeV/c² for Σ+. " \
          "M(γγ) ∈ (0.500, 0.560) GeV/c² for η(→γγ). " \
          "M(π0π0) ∉ (0.440, 0.520) GeV/c² to veto pK_S^0 background. " \
          "M_BC signal fitted with ARGUS background + MC signal shape ⊗ Gaussian.")
    .note(:reference_channel,
          "Reference decay Λc+ → Σ+ π0 uses same selection but with χ² < 17. " \
          "BF ratio = (N_sig / ε_sig) / (N_ref / ε_ref) × B_inter_ref / B_inter_sig.")
    .note(:antineutron_tag,
          "Anti-proton recoiling against detected Λc+ required to suppress combinatorial backgrounds.")
    .note(:energy_points,
          "7 energy points: 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV. " \
          "Simultaneous unbinned maximum likelihood fit to M_BC at all energy points. " \
          "ECMS constant set to 4.600 GeV as reference; per-point energy from MeasuredEcmsSvc.")

alg1.apply(event_selection_1)
alg1.execute_on(all_data + all_inc_mc + [exMC_sig_eta_gg])

# =====================================================================
# Algorithm 2: Λc+ → Σ+ η (η → π+π-π0), reference: Λc+ → Σ+ π0
# =====================================================================

decay_card_sig_eta_3pi = <<~DECAYCARD
  NoDecay e+ e-

  Decay Lambda_c+
  1.0  Sigma+ eta  PHSP;
  Enddecay

  #{sub_Sigma_p_pi0}
  #{sub_eta_3pi}
  #{sub_pi0_gg}
  End
DECAYCARD

exMC_sig_eta_3pi = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name     = "Lc_to_Sigma_eta_3pi"
  config.decay_card      = decay_card_sig_eta_3pi
  config.events          = 200_000
  config.cross_section   = :default
end

alg2 = Algorithm.new("Lc2SigmaEta3Pi")
alg2.set_header(["Lc2SigmaEta3PiAlg/Lc2SigmaEta3Pi.h"])

event_selection_2 = Selection.new
event_selection_2.select_track {
                   cos_theta 0.93
                   Vz        10.0
                   Vr        1.0
                   nChrp     ">=2"
                   nChrn     ">=2"
                 }
                .select_photon {
                   tdc_emc_start     0
                   tdc_emc_end       14
                   angle_to_track    10.0
                   energyThreshold_b 0.025
                   energyThreshold_e 0.050
                   nGam              ">=4"    # 2 from π0(Σ+) + 2 from π0(η)
                 }
                .pid(method: :probability) {
                   prob_cut 0.0
                   identify :prp, against: [:kp, :pip]
                   identify :kp,  against: :pip
                   identify :pip,  against: :kp
                 }
                .kalman_kinematic_fit([:gamma, :gamma]) {
                   invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                   chi2_cut 100
                 }
                .kinematic_fit([:prp, :pi0, :pip, :pim, :pi0]) {
                   nominal
                   constrain_four_momentum
                   invariant_mass_of(:prp, :pi0).constrain_to_nominal_mass_of(:"Sigma+")
                   chi2_cut 17
                 }

alg2.with_decay_card(decay_card_sig_eta_3pi)
    .note(:eta_3pi_selection,
          "η → π+π-π0: M(π+π-π0) ∈ (0.535, 0.560) GeV/c². " \
          "3π invariant mass constrained to η mass in ROOT fit; not expressible as a " \
          "single DSL constraint because η is a 3-body composite and the DSL's " \
          "kinematic_fit participant list does not support sub-combination mass constraints " \
          "across independent intermediate particles.")

alg2.apply(event_selection_2)
alg2.execute_on(all_data + all_inc_mc + [exMC_sig_eta_3pi])

# =====================================================================
# Algorithm 3: Λc+ → Σ+ η' (η' → π+π-η, η→γγ), reference: Λc+ → Σ+ ω
# =====================================================================

decay_card_sig_etap = <<~DECAYCARD
  NoDecay e+ e-

  Decay Lambda_c+
  1.0  Sigma+ eta'  PHSP;
  Enddecay

  Decay eta'
  1.0  pi+ pi- eta  PHSP;
  Enddecay

  #{sub_Sigma_p_pi0}
  #{sub_eta_gg}
  #{sub_pi0_gg}
  End
DECAYCARD

exMC_sig_etap = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name     = "Lc_to_Sigma_etap"
  config.decay_card      = decay_card_sig_etap
  config.events          = 200_000
  config.cross_section   = :default
end

alg3 = Algorithm.new("Lc2SigmaEtap")
alg3.set_header(["Lc2SigmaEtapAlg/Lc2SigmaEtap.h"])

event_selection_3 = Selection.new
event_selection_3.select_track {
                   cos_theta 0.93
                   Vz        10.0
                   Vr        1.0
                   nChrp     ">=2"
                   nChrn     ">=2"
                 }
                .select_photon {
                   tdc_emc_start     0
                   tdc_emc_end       14
                   angle_to_track    10.0
                   energyThreshold_b 0.025
                   energyThreshold_e 0.050
                   nGam              ">=4"
                 }
                .pid(method: :probability) {
                   prob_cut 0.0
                   identify :prp, against: [:kp, :pip]
                   identify :kp,  against: :pip
                   identify :pip,  against: :kp
                 }
                .kalman_kinematic_fit([:gamma, :gamma]) {
                   invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                   chi2_cut 100
                 }
                .kalman_kinematic_fit([:gamma, :gamma]) {
                   invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                   chi2_cut 100
                 }
                .kinematic_fit([:prp, :pi0, :pip, :pim, :eta]) {
                   nominal
                   constrain_four_momentum
                   invariant_mass_of(:prp, :pi0).constrain_to_nominal_mass_of(:"Sigma+")
                   chi2_cut 30     # optimized by FOM for Σ+η'
                 }

alg3.with_decay_card(decay_card_sig_etap)
    .note(:etap_selection,
          "η' → π+π-η: M(π+π-η) ∈ (0.946, 0.968) GeV/c². " \
          "Best candidate chosen by minimum |M(π+π-η) - M(η')|. " \
          "Reference ω → π+π-π0: M(π+π-π0) ∈ (0.760, 0.800) GeV/c², " \
          "best candidate chosen by minimum |M(π+π-π0) - M(ω)|.")

alg3.apply(event_selection_3)
alg3.execute_on(all_data + all_inc_mc + [exMC_sig_etap])