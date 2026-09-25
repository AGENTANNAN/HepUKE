# 2209.14744v2: Observation of e+e- → π0π0ψ2(3823)
# Partial reconstruction with missing photon
# Energy range: 4.23-4.70 GeV

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# === 20 energy points from Table 1 ===
data_points = [
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("705_4290"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("705_4380"),
  DatasetManager.real_data.find("705_4400"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("705_4440"),
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
]

# Decay: e+e- → π0π0ψ2(3823), ψ2(3823) → γχc1, χc1 → γJ/ψ, J/ψ → ℓ+ℓ-, π0 → γγ
decay_card = <<~DECAY
Decay psi(4260)
1.000 pi0 pi0 psi2_3823 PHSP;
Enddecay
Decay psi2_3823
1.000 gamma chi_c1 PHSP;
Enddecay
Decay chi_c1
1.000 gamma J/psi VVP;
Enddecay
Decay J/psi
1.000 e+ e- VLL;
Enddecay
Decay pi0
1.000 gamma gamma PHSP;
Enddecay
DECAY

algorithm = Algorithm.new("Pi0Pi0Psi23823", '00-00-01')
algorithm.set_header(["Pi0Pi0Psi23823Alg/Pi0Pi0Psi23823.h"])
         .set_constant({ "ECMS" => [:double, 4.360] })
         .with_decay_card(decay_card)

# Event selection: partial reconstruction with missing photon
# Require at least 5 photons (Nγ≥5), allowing one missing
# Nγ≤6 applied to suppress π0π0ψ(2S) background
event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=5"
  end
  .pid do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.1
    nlp "==1"; nlm "==1"
  end
  # 4C kinematic fit with:
  # - J/ψ mass constraint on ℓ+ℓ-
  # - Missing photon mass constrained to 0 (via miss_track_of)
  # - Two π0 mass constraints from γγ pairs
  # - χ² < 15
  .kinematic_fit([:lp, :lm, :gamma, :gamma, :gamma, :gamma, :gamma]) do
    miss_track_of :gamma
    constrain_four_momentum
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 15
    nominal
  end

algorithm.apply(event_selection)

# Note: post-kinematic-fit selections (χc1 mass window, ηJ/ψ veto,
# ψ(2S) veto, M(γγπ0π0) > 0.70 GeV/c²) are applied in the ROOT stage.
# The Nγ≤6 cut to suppress π0π0ψ(2S)→π0π0π0π0J/ψ background
# cannot be expressed via nGam alone; applied in ROOT.
algorithm.note(:post_fit_cuts, "χc1 mass window 3.49-3.53 GeV/c²; " \
  "ψ(2S) veto 3.665-3.700 GeV/c²; ηJ/ψ veto M(γγπ0π0)>0.70 GeV/c²; " \
  "Nγ≤6 suppression for π0π0π0π0J/ψ background")
algorithm.note(:chi_c1_reco, "Two remaining photons after π0 reco are boosted to " \
  "ψ2(3823) CM frame; lower-energy γ from ψ2(3823) decay, higher-energy γ + " \
  "J/ψ → χc1. ψ2(3823) mass measured via M(γγJ/ψ) ≡ M(γγℓ+ℓ-) − M(ℓ+ℓ-) + m(J/ψ)")

# Exclusive MC for signal
ex_mcs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_pi0pi0_psi23823"
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

algorithm.execute_on(data_points + ex_mcs)