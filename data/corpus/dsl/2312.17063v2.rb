# Paper 4: 2312.17063v2 — Search for a massless particle beyond the SM in Σ+ → p + invisible
# J/ψ → Σ+ anti-Σ-, Σ+ → p + invisible, anti-Σ- → anti-p π0 (ST/DT technique)
# Single energy point: √s = 3.097 GeV (J/ψ)

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card: J/ψ → Σ+ anti-Σ-
# Signal: Σ+ → p + invisible (massless particle), anti-Σ- → anti-p π0
# Note: This is a hyperon ST/DT analysis, NOT a D-tag TagAnalysis.
# The ST/DT technique is implemented via ordinary Selection with partial reconstruction.

decay_card_sigma_invisible = <<~DECAYCARD
  Decay J/psi
  1 Sigma+ anti-Sigma- PHSP;
  Enddecay
  Decay Sigma+
  1 p+ invisible PHSP;
  Enddecay
  Decay anti-Sigma-
  1 anti-p- pi0 PHSP;
  Enddecay
  Decay pi0
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

sig_sigma_inv = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_sigma_invisible"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_sigma_invisible
  config.cross_section   = :default
end

# ========== Algorithm ==========
alg_sigma_inv = Algorithm.new("SigmaInvisible", version: '00-00-01')
alg_sigma_inv.set_header(["SigmaInvisibleAlg/SigmaInvisible.h"])
alg_sigma_inv.set_constant({ "ECMS" => [:double, 3.097] })

# ========== ST Selection: anti-Σ- → anti-p π0 ==========
# This is an ordinary reconstruction using Selection, not TagAnalysis.
# The anti-Σ- is reconstructed from its decay products: anti-proton + π0(→γγ)
# π0 → γγ, constrained via Kalman fit.
# The ST yield is obtained from the M_BC distribution

st_selection = Selection.new
  # Track selection
  .select_track do
    nTot ">=2"
    cos_theta 0.93
    Vz 10.0
    Vr 2.0
  end
  # Photon selection: >=2 photons for π0 reconstruction
  .select_photon do
    nGam ">=2"
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
  end
  # Photon selection: isolation from anti-proton (20° for anti-proton nuclear interactions)
  .select_isolated_photon do
    angle_to_prm_track 20.0
    angle_to_prp_track 10.0
  end
  # PID: identify protons
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pip, :kp]
  end
  # π0 reconstruction from γγ pairs (1C Kalman fit, χ²<25)
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
  end
  # Build anti-Σ- from anti-p π0
  .build_virtual_particle(:anti_Sigma_minus, from: [:proton_m, :pi0])

# Note: ROOT-level:
# - π0 mass window: (115, 150) MeV/c²
# - Best anti-Σ- combination: closest M(anti-p π0) to nominal anti-Σ- mass
# - M(anti-p π0) window: |M - M_Σ| < 15 MeV/c²
# - M_BC fit for ST yield extraction

# DT event selection is complex and uses:
# - Exactly 1 additional proton
# - 2C kinematic fit: J/ψ → p anti-p π0 + invisible (constrain π0 mass, invisible mass=0)
# - Competing 2C fit: invisible mass constrained to π0 mass (veto Σ+→pπ0)
# - 5C kinematic fit: additional photon on DT side (veto Σ+→pγ)
# - 6C kinematic fit: additional π0 on DT side (veto Σ+→pπ0)
# - Primary and secondary vertex fits, decay length L/σL > 2
# - cos θ_invisible cut: |cos θ_inv| < 0.8
# - M(p+inv) in [1.18, 1.20] GeV/c²
# These are all ROOT-level operations due to their complexity.

alg_sigma_inv.with_decay_card(decay_card_sigma_invisible).apply(st_selection)
alg_sigma_inv.execute_on([jpsi_data, jpsi_incMC, sig_sigma_inv])

# Note: This analysis uses a hyperon ST/DT technique, which differs from the D-tag/Λc-tag
# TagAnalysis framework. The ST reconstruction and DT selection are implemented via ordinary
# Selection with partial reconstruction and ROOT-level kinematic fits.
# Key ROOT-level procedures:
# - ST M_BC distribution fit to extract ST yield
# - DT: 2C, 5C, 6C kinematic fits with competing hypotheses
# - Vertex fits: primary vertex + secondary vertex for p anti-p
# - Data-driven E_extra correction for anti-proton nuclear interactions
# - Binned maximum-likelihood fit to E_extra for DT signal yield
# - Bayesian upper limit on BF(Σ+ → p + invisible)