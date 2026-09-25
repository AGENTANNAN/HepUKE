# BESIII ψ(3686) → φ η η' PWA analysis
# ArXiv: 2410.05736v1  —  Observation of X(2300) axial-vector state
# Dataset: (2712.4 ± 14.3)×10^6 ψ(3686) events

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: ψ(3686) → φ η η'
# φ → K+K-, η → γγ, η' → γ π+π-
decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000  phi eta etap  PHSP;
    Enddecay

    Decay phi
    1.0000  K+  K-  VSS;
    Enddecay

    Decay eta
    1.0000  gamma  gamma  PHSP;
    Enddecay

    Decay etap
    1.0000  gamma  pi+  pi-  PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_phi_eta_etap"
  config.related_dataset = psip_data
  config.events = 1_000_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Algorithm ###
alg = Algorithm.new("PsipPhiEtaEtap")
alg.set_header(["PsipPhiEtaEtapAlg/PsipPhiEtaEtap.h"])
  .set_constant(ECMS: 3.686)

### Event selection (BOSS) — stops at the 5C kinematic fit ###
event_selection = Selection.new
  .select_track {
    cos_theta  0.93
    Vz         10.0
    Vr         1.0
    nChrp      "==2"
    nChrn      "==2"
    nNet       "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=3"
  }
  .pid(method: :probability) {
    prob_cut   0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  }
  # 4C kinematic fit (pre-selection): ψ(3686) → K+K- π+π- γγγ
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
    chi2_cut 50
  }
  # Competing-hypothesis veto: 4-photon background (ψ(3686) → K+K- π+π- γγγγ)
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }
  # Competing-hypothesis veto: 2-photon background (ψ(3686) → K+K- π+π- γγ)
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {
    constrain_four_momentum
  }
  # 5C kinematic fit (nominal): 4C + η mass constraint
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  }

alg.note(:background_veto, "chi_cJ veto: |M(K+K-pi+pi-gamma gamma) - M(chi_cJ)| > 0.010 GeV/c^2; J/psi veto: |M(gamma gamma gamma K+K-) - M(J/psi)| > 0.030 GeV/c^2; eta J/psi veto: |M(gamma pi+pi-K+K-) - M(J/psi)| > 0.030 GeV/c^2")
  .note(:helix_correction, "helix parameter correction applied to simulated charged tracks to match data resolution in the 4C kinematic fit")
  .with_decay_card(decay_card)
  .apply(event_selection)

alg.execute_on([psip_data, psip_incMC, exMC_signal])