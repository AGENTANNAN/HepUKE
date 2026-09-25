# 2210.09601v2: First observation of J/ψ/ψ(3686) → ηΣ+Σ-
# Two separate datasets: J/ψ (3.097 GeV) and ψ(3686) (3.686 GeV)
# Decay chain: η → γγ, Σ+ → pπ0, Σ- → pbar π0
# 7C kinematic fit: 4C + η mass + 2×π0 mass

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_jpsi    = DatasetManager.real_data.find("708_3097")
data_psip    = DatasetManager.real_data.find("709_3686")
inc_mc_jpsi  = DatasetManager.inclusive_mc.find("708_3097")
inc_mc_psip  = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: ψ → η Σ+ Σ-, with η → γγ, Σ+ → pπ0, Σ- → pbar π0, π0 → γγ
decay_card_jpsi = <<~DECAY
Decay J/psi
1.000 eta Sigma+ anti-Sigma- PHSP;
Enddecay
Decay eta
1.000 gamma gamma PHSP;
Enddecay
Decay Sigma+
1.000 p+ pi0 PHSP;
Enddecay
Decay anti-Sigma-
1.000 anti-p- pi0 PHSP;
Enddecay
Decay pi0
1.000 gamma gamma PHSP;
Enddecay
DECAY

decay_card_psip = <<~DECAY
Decay psi(2S)
1.000 eta Sigma+ anti-Sigma- PHSP;
Enddecay
Decay eta
1.000 gamma gamma PHSP;
Enddecay
Decay Sigma+
1.000 p+ pi0 PHSP;
Enddecay
Decay anti-Sigma-
1.000 anti-p- pi0 PHSP;
Enddecay
Decay pi0
1.000 gamma gamma PHSP;
Enddecay
DECAY

# ======================================================================
# J/ψ → ηΣ+Σ-
# 7C fit: χ²_7C < 30;  η' veto applied
# Σ+Σ- pairing via minimizing Δ = √((Mpπ0−mΣ+)² + (Mpbarπ0−mΣ-)²)
# ======================================================================
algorithm_jpsi = Algorithm.new("EtaSigmaSigmaJpsi", '00-00-01')
algorithm_jpsi.set_header(["EtaSigmaSigmaJpsiAlg/EtaSigmaSigmaJpsi.h"])
               .set_constant({ "ECMS" => [:double, 3.097] })
               .with_decay_card(decay_card_jpsi)

selection_jpsi = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        20.0
    nChrp     "==1"
    nChrn     "==1"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=6"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
    nprp "==1"; nprm "==1"
  end
  # 7C kinematic fit:
  #   4C (momentum conservation)
  #   + η mass constraint (+1C)
  #   + π0 mass constraint (+1C)
  #   + π0 mass constraint (+1C)
  #   = 7C total
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) do
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 30
    nominal
  end

algorithm_jpsi.apply(selection_jpsi)

algorithm_jpsi.note(:post_fit_cuts,
  "η mass window: (0.450, 0.650) GeV/c². " \
  "π0 mass window: (0.110, 0.160) GeV/c². " \
  "Σ+ signal region: M(pπ0) ∈ [1.177, 1.201] GeV/c². " \
  "Σ- signal region: M(pbarπ0) ∈ [1.177, 1.201] GeV/c². " \
  "Σ sidebands: [1.141, 1.165] ∪ [1.213, 1.237] GeV/c². " \
  "Σ+Σ- pairing: minimize Δ = √((Mpπ0−mΣ+)² + (Mpbarπ0−mΣ-)²).")

algorithm_jpsi.note(:eta_prime_veto,
  "η' → ηπ0π0 background veto: require " \
  "M(ηπ0π0) ∉ [0.95, 0.97] GeV/c² to suppress J/ψ → p pbar η'.")

algorithm_jpsi.note(:peaking_background,
  "Peaking background from J/ψ → π0Σ+Σ- (107.6±0.6 events) " \
  "subtracted using control sample and exclusive MC shapes.")

# ======================================================================
# ψ(3686) → ηΣ+Σ-
# 7C fit: χ²_7C < 25;  recoil mass cut applied
# ======================================================================
algorithm_psip = Algorithm.new("EtaSigmaSigmaPsip", '00-00-01')
algorithm_psip.set_header(["EtaSigmaSigmaPsipAlg/EtaSigmaSigmaPsip.h"])
               .set_constant({ "ECMS" => [:double, 3.686] })
               .with_decay_card(decay_card_psip)

selection_psip = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        20.0
    nChrp     "==1"
    nChrn     "==1"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=6"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
    nprp "==1"; nprm "==1"
  end
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) do
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    nominal
  end

algorithm_psip.apply(selection_psip)

algorithm_psip.note(:post_fit_cuts,
  "η mass window: (0.450, 0.650) GeV/c². " \
  "π0 mass window: (0.110, 0.160) GeV/c². " \
  "Σ+ signal region: M(pπ0) ∈ [1.177, 1.201] GeV/c². " \
  "Σ- signal region: M(pbarπ0) ∈ [1.177, 1.201] GeV/c². " \
  "Σ sidebands: [1.141, 1.165] ∪ [1.213, 1.237] GeV/c².")

algorithm_psip.note(:recoil_mass_cuts,
  "ψ(3686) → ηJ/ψ and ψ(3686) → γχc0,1,2 background suppressed via: " \
  "Mrec(η) < 3.050 GeV/c² (recoil mass of η). " \
  "ψ(3686) → π0π0J/ψ veto: Mrec(π0π0) ∉ [3.080, 3.120] GeV/c².")

algorithm_psip.note(:peaking_background,
  "Peaking background from ψ(3686) → γχc0,1,2, χc0,1,2 → π0Σ+Σ- " \
  "(1.1±0.1 events) subtracted using control sample and exclusive MC shapes.")

# Exclusive MC for J/ψ signal
ex_mc_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_eta_sigma_sigma_jpsi"
  config.related_dataset = data_jpsi
  config.events          = 500_000
  config.decay_card      = decay_card_jpsi
  config.cross_section   = :default
end

# Exclusive MC for ψ(3686) signal
ex_mc_psip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_eta_sigma_sigma_psip"
  config.related_dataset = data_psip
  config.events          = 500_000
  config.decay_card      = decay_card_psip
  config.cross_section   = :default
end

algorithm_jpsi.execute_on([data_jpsi, inc_mc_jpsi, ex_mc_jpsi])
algorithm_psip.execute_on([data_psip, inc_mc_psip, ex_mc_psip])