# Observation of hc radiative decays to multiple light hadrons via ψ(3686)→π0 hc
# ψ(3686), 2.712×10^8 events, arXiv:2501.15447v2
# 4 decay modes: I: γ π+π-, II: γ π+π-η, III: γ 2(π+π-), IV: γ pp̄
# Rule T1: each independent decay mode has its own Algorithm

psi3686_data = DatasetManager.real_data.find("709_3686")
psi3686_incMC = DatasetManager.inclusive_mc.find("709_3686")

### ── Mode I: hc → γ π+π- ── ###

decay_card_I = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi0 hc HELAMP 1 0 1 0 1 0 1 0;
  Enddecay

  Decay hc
  1.000 gamma pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_I = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi2s_pi0_hc_gamma_pipi"
  config.related_dataset = psi3686_data
  config.events = 100_000
  config.decay_card = decay_card_I
  config.cross_section = :default
end

alg_I = Algorithm.new("Psi2SPi0HcGammaPiPi")
alg_I.set_header(["Psi2SPi0HcGammaPiPiAlg/Psi2SPi0HcGammaPiPi.h"])
     .set_constant({ "ECMS" => [:double, 3.686] })

sel_I = Selection.new
sel_I
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrg "==2"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=3"
  }
  .pid(method: :probability) {
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  # Reconstruct π0 → γγ via 1C Kalman fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # 5C kinematic fit: 4C + π0 mass constraint
  .kinematic_fit([:pip, :pim, :gamma, :pi0]) do
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  end

alg_I.with_decay_card(decay_card_I).apply(sel_I)

alg_I.note(:root_cuts, "ROOT-level: π0 mass [0.10,0.16] GeV/c², 5C χ² min over photon combinations, cosθ(π+π-)>-0.98, ω/φ veto on M(π+π-π0) by 3σ, ω veto on M(π0γ), π0 veto on M(γ_Hγ), M(π+π-)<2.8 GeV/c², M(X)<2.8 GeV/c² for η_c veto, photon multiplicity veto χ²_4C")
alg_I.note(:signal_extraction, "Signal yield from fit to M(γπ+π-); ARGUS background; f2(1270) simultaneous fit")

### ── Mode II: hc → γ π+π-η (η→γγ) ── ###

decay_card_II = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi0 hc HELAMP 1 0 1 0 1 0 1 0;
  Enddecay

  Decay hc
  1.000 gamma pi+ pi- eta PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_II = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi2s_pi0_hc_gamma_pipi_eta"
  config.related_dataset = psi3686_data
  config.events = 100_000
  config.decay_card = decay_card_II
  config.cross_section = :default
end

alg_II = Algorithm.new("Psi2SPi0HcGammaPiPiEta")
alg_II.set_header(["Psi2SPi0HcGammaPiPiEtaAlg/Psi2SPi0HcGammaPiPiEta.h"])
      .set_constant({ "ECMS" => [:double, 3.686] })

sel_II = Selection.new
sel_II
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrg "==2"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=5"
  }
  .pid(method: :probability) {
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  # Reconstruct π0 → γγ via 1C Kalman fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # Reconstruct η → γγ via 1C Kalman fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  # 5C kinematic fit: 4C + π0 mass constraint
  .kinematic_fit([:pip, :pim, :gamma, :pi0, :eta]) do
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  end

alg_II.with_decay_card(decay_card_II).apply(sel_II)

alg_II.note(:root_cuts, "ROOT-level: π0 mass [0.10,0.16] GeV/c², η mass [0.513,0.583] GeV/c², 5C χ² min, photon multiplicity veto χ²_4C, M(π+π-η)>1.0 GeV/c² (η' veto), M(X)<2.8 GeV/c² (η_c veto), M(γγ)_mix outside π0 mass, M_recoil(π+π-)<3.085 GeV/c²")

### ── Mode III: hc → γ 2(π+π-) ── ###

decay_card_III = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi0 hc HELAMP 1 0 1 0 1 0 1 0;
  Enddecay

  Decay hc
  1.000 gamma pi+ pi- pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_III = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi2s_pi0_hc_gamma_4pi"
  config.related_dataset = psi3686_data
  config.events = 100_000
  config.decay_card = decay_card_III
  config.cross_section = :default
end

alg_III = Algorithm.new("Psi2SPi0HcGamma4Pi")
alg_III.set_header(["Psi2SPi0HcGamma4PiAlg/Psi2SPi0HcGamma4Pi.h"])
       .set_constant({ "ECMS" => [:double, 3.686] })

sel_III = Selection.new
sel_III
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrg "==4"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=3"
  }
  .pid(method: :probability) {
    identify :pion, against: [:kaon]
    npip "==2"
    npim "==2"
  }
  # Reconstruct π0 → γγ via 1C Kalman fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # 5C kinematic fit: 4C + π0 mass constraint
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :pi0]) do
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  end

alg_III.with_decay_card(decay_card_III).apply(sel_III)

alg_III.note(:root_cuts, "ROOT-level: π0 mass [0.10,0.16] GeV/c², 5C χ² min, photon multiplicity veto χ²_4C, M(π+π-π0) outside η/ω mass regions, M_recoil(π+π-) outside J/ψ [3.091,3.103] GeV/c², M(X)<2.8 GeV/c² (η_c veto)")
alg_III.note(:intermediate, "Dominated by hc→γρ0ρ0 cascade decay; simulated accordingly for efficiency")

### ── Mode IV: hc → γ pp̄ ── ###

decay_card_IV = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi0 hc HELAMP 1 0 1 0 1 0 1 0;
  Enddecay

  Decay hc
  1.000 gamma p+ anti-p- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_IV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi2s_pi0_hc_gamma_ppbar"
  config.related_dataset = psi3686_data
  config.events = 100_000
  config.decay_card = decay_card_IV
  config.cross_section = :default
end

alg_IV = Algorithm.new("Psi2SPi0HcGammaPPbar")
alg_IV.set_header(["Psi2SPi0HcGammaPPbarAlg/Psi2SPi0HcGammaPPbar.h"])
      .set_constant({ "ECMS" => [:double, 3.686] })

sel_IV = Selection.new
sel_IV
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrg "==2"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=3"
  }
  # PID: proton likelihood > kaon and pion likelihoods
  .pid(method: :probability) {
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  }
  # Reconstruct π0 → γγ via 1C Kalman fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # 5C kinematic fit: 4C + π0 mass constraint
  .kinematic_fit([:prp, :prm, :gamma, :pi0]) do
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  end

alg_IV.with_decay_card(decay_card_IV).apply(sel_IV)

alg_IV.note(:root_cuts, "ROOT-level: π0 mass [0.10,0.16] GeV/c², 5C χ² min, PID L(p)>L(K) & L(p)>L(π), photon multiplicity veto χ²_4C, M(X)<2.8 GeV/c² (η_c veto)")
alg_IV.note(:signal_extraction, "Signal yields from fits to M(γX) distributions; ARGUS smooth background; peaking backgrounds fixed to MC estimates; significances from Δ(log L) with Δndf")

alg_I.execute_on([psi3686_data, psi3686_incMC, exMC_I])
alg_II.execute_on([psi3686_data, psi3686_incMC, exMC_II])
alg_III.execute_on([psi3686_data, psi3686_incMC, exMC_III])
alg_IV.execute_on([psi3686_data, psi3686_incMC, exMC_IV])