# =====================================================================
# BESIII: ψ(3686) → γ χ_cJ, χ_cJ → 3 K_S^0 K^± π^∓
#   [arXiv:2505.15620v1]
#
# First observation of χ_c0,1,2 → 3 K_S^0 K^± π^∓ using
# (2712.4 ± 14.3) × 10^6 ψ(3686) events.
# Ordinary analysis (Algorithm + Selection).
# =====================================================================

psi2s_data  = DatasetManager.real_data.find("709_3686")
psi2s_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
  Decay psi(2S)
  0.09  gamma chi_c0  VSP_PWAVE;
  0.09  gamma chi_c1  VSP_PWAVE;
  0.09  gamma chi_c2  VSP_PWAVE;
  Enddecay

  Decay chi_c0
  1.0  K_S0 K_S0 K_S0 K+ pi-  PHSP;
  Enddecay

  Decay chi_c1
  1.0  K_S0 K_S0 K_S0 K+ pi-  PHSP;
  Enddecay

  Decay chi_c2
  1.0  K_S0 K_S0 K_S0 K+ pi-  PHSP;
  Enddecay

  Decay K_S0
  1.0  pi+ pi-  PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "chicJ_to_3KsKPi"
  config.related_dataset = psi2s_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

alg_name = "ChicJ2ThreeKsKPi"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz        10.0
                  Vr        1.0
                  nTot      ">=8"      # ≥8 charged tracks
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam              ">=1"    # ≥1 photon from ψ(3686) radiative decay
                }
               .pid(method: :probability) {
                  prob_cut 0.0
                  identify :kp, against: :pim
                  identify :km, against: :pip
                  identify :pip, against: :kp
                  identify :pim, against: :km
                }
               .secondary_vertex_fit(:K_S0, daughters: [:pip, :pim]) {
                  chi2_cut 200              # vertex fit quality < 200
                }
               .kinematic_fit([:gamma, :K_S0, :K_S0, :K_S0, :kp, :pim]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 50               # optimized by FOM
                }

algorithm
  .note(:ks_selection,
        "K_S^0 candidates reconstructed from two oppositely charged tracks with |Vz|<20 cm, " \
        "assigned as π+π- without PID. Common vertex fit χ²<200, decay length > 2σ from IP. " \
        "Signal region: |M(π+π-) - 0.498| < 0.012 GeV/c² (3σ). " \
        "Sideband regions: 0.020 < |M-0.498| < 0.044 GeV/c² (5σ-11σ) used for " \
        "background subtraction in ROOT.")
  .note(:sideband_method,
        "K_S^0 sideband method: N_net = N_sig - 1/2 N_SB1 + 1/4 N_SB2 - 1/8 N_SB3. " \
        "Simultaneous unbinned ML fit to M(3K_S^0 K± π∓) distributions in signal and " \
        "sideband regions. Signal shape: Breit-Wigner ⊗ Gaussian. " \
        "Background: 1st-order Chebyshev polynomial.")
  .note(:mixed_mc,
        "Mixed MC with sub-resonance fractions from data: " \
        "2K_S^0 K*(892)K (33-42%), 3K_S^0 K*(892)^0 + c.c. (22-31%), " \
        "f_X(1500) K_S^0 Kπ (21-33%), and non-resonant PHSP (3-10%), " \
        "as determined in Table II of the paper.")
  .note(:radiative_photon,
        "The radiative photon from ψ(3686) → γ χ_cJ is the highest-energy photon candidate. " \
        "Best combination chosen by minimum χ²_4C.")
  .note(:process,
        "Decay chain: ψ(3686) → γ χ_cJ, χ_cJ → 3K_S^0 K± π∓, K_S^0 → π+π-. " \
        "χ²_4C < 50 chosen by FOM optimization. The 3 K_S^0 candidates must be in the " \
        "signal region for N_sig, with sideband regions for combinatorial background subtraction.")

algorithm.with_decay_card(decay_card_signal).apply(event_selection)
algorithm.execute_on([psi2s_data, psi2s_incMC, exMC_signal])