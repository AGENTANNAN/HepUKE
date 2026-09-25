# BOSS Ruby DSL for BESIII paper 1802.00583v3
# Observation of a0(980)-f0(980) mixing
# Two main analyses:
#   Analysis A: J/ψ → φ η π0 (φ → K+K-, η → γγ or π+π-π0)
#   Analysis B: ψ(3686) → γ χc1, χc1 → π0 π+π- (mixing signal)

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card A1: J/ψ → φ η π0, φ → K+K-, η → γγ, π0 → γγ
decay_card_A1 = <<~DECAYCARD
    Decay J/psi
    1.000  phi  eta  pi0  PHSP;
    Enddecay

    Decay phi
    1.000  K+  K-  VSS;
    Enddecay

    Decay eta
    1.000  gamma  gamma  PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card A2: J/ψ → φ η π0, φ → K+K-, η → π+π-π0, π0 → γγ
decay_card_A2 = <<~DECAYCARD
    Decay J/psi
    1.000  phi  eta  pi0  PHSP;
    Enddecay

    Decay phi
    1.000  K+  K-  VSS;
    Enddecay

    Decay eta
    1.000  pi+  pi-  pi0  PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card B: ψ(3686) → γ χc1, χc1 → π0 π+π-
decay_card_B = <<~DECAYCARD
    Decay psi(2S)
    1.000  gamma  chi_c1  PHSP;
    Enddecay

    Decay chi_c1
    1.000  pi+  pi-  pi0  PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples
exMC_A1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_jpsi_phi_eta_pi0_eta2gg"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_A1
  config.cross_section = :default
end

exMC_A2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_jpsi_phi_eta_pi0_eta2pipipi0"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_A2
  config.cross_section = :default
end

exMC_B = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_psip_gamma_chic1_pipipi0"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_B
  config.cross_section = :default
end

### Event selection - Mode A1: J/ψ → φ η π0 (η → γγ, φ → K+K-, π0 → γγ) ###
# Required: 2 kaons (opposite charge), >= 4 photons
# 4C kinematic fit: K+K-γγγγ, χ² < 50

alg_A1 = Algorithm.new("JpsiPhiEtaPi0_eta2gg")
alg_A1.set_header(["JpsiPhiEtaPi0_eta2ggAlg/JpsiPhiEtaPi0_eta2gg.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

sel_A1 = Selection.new
sel_A1.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=1"
    nNet  "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=4"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"; nkm ">=1"
  end
  .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma]) do
    nominal
    constrain_four_momentum
    chi2_cut 50
  end

alg_A1
  .note(:pid_correction_method, "tracking and PID efficiency corrections: 1% per track")
  .note(:helix_correction, "kinematic fit efficiency correction: 1.5% from control sample")
  .note(:efficiency_curve, "post-fit π0 and η selection via χ² minimization: |M(γγ)-m_π0| < 15 MeV/c², |M(γγ)-m_η| < 30 MeV/c²; φ mass window |M(K+K-)-m_φ| < 10 MeV/c²; sideband subtraction for backgrounds")
  .note(:background_veto, "veto J/ψ→K+K-π0π0 via χ²(π0π0) > 40 and J/ψ→K+K-ηη via χ²(ηη) > 5; applied at post-fit level")
  .with_decay_card(decay_card_A1)
  .apply(sel_A1)

### Event selection - Mode A2: J/ψ → φ η π0 (η → π+π-π0, φ → K+K-) ###
# Required: 2 kaons + 2 pions (opposite charge each), >= 4 photons
# 4C kinematic fit: K+K-π+π-γγγγ, χ² < 60

alg_A2 = Algorithm.new("JpsiPhiEtaPi0_eta2pipipi0")
alg_A2.set_header(["JpsiPhiEtaPi0_eta2pipipi0Alg/JpsiPhiEtaPi0_eta2pipipi0.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

sel_A2 = Selection.new
sel_A2.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
    nNet  "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=4"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"; nkm ">=1"
    npip ">=1"; npim ">=1"
  end
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma, :gamma]) do
    nominal
    constrain_four_momentum
    chi2_cut 60
  end

alg_A2
  .note(:pid_correction_method, "tracking and PID efficiency corrections: 1% per track")
  .note(:helix_correction, "kinematic fit efficiency correction: 1.5% from control sample")
  .note(:efficiency_curve, "post-fit η, π0, φ selection via mass windows: η mass window 20 MeV/c², π0 mass window 15 MeV/c², φ mass window 10 MeV/c²")
  .note(:background_veto, "φ and η sidebands used for background estimation")
  .with_decay_card(decay_card_A2)
  .apply(sel_A2)

### Event selection - Mode B: ψ(3686) → γ χc1, χc1 → π0 π+π- ###
# Required: 2 pions (opposite charge), >= 3 photons
# 4C kinematic fit: π+π-γγγ, χ² < 20
# Veto: χ²_4C(π+π-γγγ) < χ²_4C(π+π-γγγγ) and χ²_4C(π+π-γγγ) < χ²_4C(π+π-γγ)

alg_B = Algorithm.new("PsipChic1Pi0PiPi")
alg_B.set_header(["PsipChic1Pi0PiPiAlg/PsipChic1Pi0PiPi.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

sel_B = Selection.new
sel_B.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=1"
    nNet  "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=3"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"; npim ">=1"
  end
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma]) do
    nominal
    constrain_four_momentum
    chi2_cut 20
  end
  # Competing-hypothesis veto: π+π-4γ hypothesis (no chi2_cut, no nominal)
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma]) do
    constrain_four_momentum
  end
  # Competing-hypothesis veto: π+π-2γ hypothesis (no chi2_cut, no nominal)
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) do
    constrain_four_momentum
  end

alg_B
  .note(:pid_correction_method, "tracking and PID efficiency corrections: 1% per track; photon detection 1% per photon")
  .note(:helix_correction, "kinematic fit efficiency correction: 2.5% from control sample ψ(3686)→π+π-J/ψ→π+π-γη")
  .note(:efficiency_curve, "post-fit χc1 mass window: |M(π+π-π0)-m_χc1| < 20 MeV/c²; π0 mass window 15 MeV/c²; competing-hypothesis vetoes via chi2 comparison")
  .with_decay_card(decay_card_B)
  .apply(sel_B)

### Execute ###
alg_A1.execute_on([jpsi_data, jpsi_incMC, exMC_A1])
alg_A2.execute_on([jpsi_data, jpsi_incMC, exMC_A2])
alg_B.execute_on([psip_data, psip_incMC, exMC_B])