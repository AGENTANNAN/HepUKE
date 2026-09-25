# Paper: 2203.10439v1 — Search for new hadronic decays of hc and observation of hc → pp̅η
# Process: ψ(3686) → π0 hc, hc → pp̅X (X = π+π-π0, η, π0)
# Three independent decay modes → three separate Algorithm objects (Rule T1)
# Single energy: ψ(3686) at 3.686 GeV, ORDINARY analysis

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ── Datasets ──────────────────────────────────────────────────────────
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ── Common decay card: ψ(3686) → π0 hc ──────────────────────────────
# All three modes share the same production: ψ(3686) → π0 hc
# The hc decay differs per mode

# ═══════════════════════════════════════════════════════════════════════
# Mode I: hc → pp̅ π+π-π0
# ═══════════════════════════════════════════════════════════════════════
decay_card_I = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi0 hc  PHSP;
  Enddecay
  Decay hc
  1.0000 p+ anti-p- pi+ pi- pi0  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

ex_mc_I = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_hc_pp_pipip_pi0"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_I
  config.cross_section   = :default
end

algo_I = Algorithm.new("HcToPPbarPipPimPi0")

event_sel_I = Selection.new
  .select_track do
    nChrp ">=2"
    nChrn ">=2"
    nTot  ">=4"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .select_photon do
    nGam              ">=6"    # 2π0 → 4γ + bachelor π0 → 2γ = 6γ
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
  end
  .pid(method: :probability) do
    identify(:prp, :prm, against: :pion)
    identify(:pip, :pim, against: :kaon)
    prob_cut 0.001
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=3"
  end
  .kinematic_fit([:prp, :prm, :pip, :pim, :pi0, :pi0, :pi0]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # bachelor π0 from ψ(3686) decay (if needed)
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algo_I
  .set_header(["HcToPPbarPipPimPi0/HcToPPbarPipPimPi0.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })
  .note(:helix_correction, "track helix parameter correction applied before the 6C kinematic fit; difference with/without corrections taken as systematic uncertainty")
  .note(:background_veto, "vetoes applied at ROOT level: χ2_4C,nγ < χ2_4C,(n-1)γ and χ2_4C,nγ < χ2_4C,(n+1)γ (photon number test); |RM(π+π-)−m(J/ψ)| < 15-16 MeV/c2; |RM(π0π0)−m(J/ψ)| < 16 MeV/c2; |M(π+π-π0)−m(η)| < 6 MeV/c2; ω veto M(π+π-π0) ∈ (762,802) MeV/c2; Σ± veto M(pπ0) ∈ (1180,1196) and M(p̅π0) ∈ (1181,1194) MeV/c2; Λ/Λ̅ veto M(pπ-) ∈ (1104,1122) and M(p̅π+) ∈ (1104,1122) MeV/c2; Ks0 veto M(π+π-) ∈ (490,499) MeV/c2")
  .note(:"6C_kinematic_fit", "6C kinematic fit: four-momentum conservation (4C) + two π0 mass constraints (2C); the bachelor π0 from ψ(3686)→π0hc is included as one of the three π0 candidates; χ2_6C cut at 45 optimized via FOM")
  .note(:bachelor_pi0_identification, "the bachelor π0 from ψ(3686)→π0hc is identified as the π0 whose energy is closest to the expected energy in the ψ(3686) rest frame")
  .with_decay_card(decay_card_I)
  .apply(event_sel_I)

algo_I.execute_on([psip_data, psip_incMC, ex_mc_I])

# ═══════════════════════════════════════════════════════════════════════
# Mode II: hc → pp̅ η  (with η → γγ and η → π+π-π0 sub-modes)
# ═══════════════════════════════════════════════════════════════════════

# Mode IIa: η → γγ
decay_card_IIa = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi0 hc  PHSP;
  Enddecay
  Decay hc
  1.0000 p+ anti-p- eta  PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

ex_mc_IIa = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_hc_pp_eta_gg"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_IIa
  config.cross_section   = :default
end

algo_IIa = Algorithm.new("HcToPPbarEtaGG")

event_sel_IIa = Selection.new
  .select_track do
    nChrp ">=1"
    nChrn ">=1"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .select_photon do
    nGam              ">=4"    # bachelor π0 → 2γ + η → 2γ = 4γ
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
  end
  .pid(method: :probability) do
    identify(:prp, :prm, against: :pion)
    prob_cut 0.001
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  end
  .kinematic_fit([:prp, :prm, :pi0, :eta]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algo_IIa
  .set_header(["HcToPPbarEtaGG/HcToPPbarEtaGG.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })
  .note(:helix_correction, "track helix parameter correction applied before the kinematic fit")
  .note(:background_veto, "vetoes: |RM(η)−m(J/ψ)| < 22 MeV/c2; |RM(π0π0)−m(J/ψ)| < 30 MeV/c2; Σ± veto M(pπ0) ∈ (1172,1202) and M(p̅π0) ∈ (1171,1202) MeV/c2; photon-number test as in Mode I; η_c peaking background from hc→γη_c, η_c→pp̅π0 included in simultaneous fit")
  .note(:"6C_kinematic_fit", "6C kinematic fit: 4C + π0 mass (1C) + η mass (1C); χ2_6C < 45; the bachelor π0 from ψ(3686)→π0hc is identified by energy")
  .note(:simultaneous_fit, "simultaneous unbinned maximum likelihood fit to pp̅η spectra from both η→γγ and η→π+π-π0 sub-modes; common B(hc→pp̅η) parameter")
  .with_decay_card(decay_card_IIa)
  .apply(event_sel_IIa)

algo_IIa.execute_on([psip_data, psip_incMC, ex_mc_IIa])

# Mode IIb: η → π+π-π0
decay_card_IIb = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi0 hc  PHSP;
  Enddecay
  Decay hc
  1.0000 p+ anti-p- eta  PHSP;
  Enddecay
  Decay eta
  1.0000 pi+ pi- pi0  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

ex_mc_IIb = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_hc_pp_eta_3pi"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_IIb
  config.cross_section   = :default
end

algo_IIb = Algorithm.new("HcToPPbarEta3Pi")

event_sel_IIb = Selection.new
  .select_track do
    nChrp ">=2"
    nChrn ">=2"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .select_photon do
    nGam              ">=4"    # bachelor π0 → 2γ + η sub-mode π0 → 2γ = 4γ
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
  end
  .pid(method: :probability) do
    identify(:prp, :prm, against: :pion)
    identify(:pip, :pim, against: :kaon)
    prob_cut 0.001
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=2"    # bachelor π0 + π0 from η decay
  end
  .build_virtual_particle(:eta, from: [:pip, :pim, :pi0])
  .kinematic_fit([:prp, :prm, :pi0, :eta]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algo_IIb
  .set_header(["HcToPPbarEta3Pi/HcToPPbarEta3Pi.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })
  .note(:helix_correction, "track helix parameter correction applied before the kinematic fit")
  .note(:background_veto, "same veto scheme as Mode IIa; η signal region M(π+π-π0) ∈ (532, 562) MeV/c2")
  .note(:"6C_kinematic_fit", "6C kinematic fit: 4C + bachelor π0 mass (1C) + η from π+π-π0 mass (1C); stored as virtual particle :eta")
  .note(:simultaneous_fit, "simultaneous fit with Mode IIa for common B(hc→pp̅η)")
  .with_decay_card(decay_card_IIb)
  .apply(event_sel_IIb)

algo_IIb.execute_on([psip_data, psip_incMC, ex_mc_IIb])

# ═══════════════════════════════════════════════════════════════════════
# Mode III: hc → pp̅ π0
# ═══════════════════════════════════════════════════════════════════════
decay_card_III = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi0 hc  PHSP;
  Enddecay
  Decay hc
  1.0000 p+ anti-p- pi0  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

ex_mc_III = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_hc_pp_pi0"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_III
  config.cross_section   = :default
end

algo_III = Algorithm.new("HcToPPbarPi0")

event_sel_III = Selection.new
  .select_track do
    nChrp ">=1"
    nChrn ">=1"
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .select_photon do
    nGam              ">=4"    # 2 π0 → 4γ
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
  end
  .pid(method: :probability) do
    identify(:prp, :prm, against: :pion)
    prob_cut 0.001
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=2"
  end
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algo_III
  .set_header(["HcToPPbarPi0/HcToPPbarPi0.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })
  .note(:helix_correction, "track helix parameter correction applied before the 6C kinematic fit")
  .note(:background_veto, "vetoes: |RM(π0π0)−m(J/ψ)| < 16-30 MeV/c2; Σ± and other vetoes as in Table I for Mode III; photon-number test; χ2_6C < 64")
  .note(:"6C_kinematic_fit", "6C kinematic fit: 4C + two π0 mass constraints (2C); bachelor π0 identified by energy proximity to expected value")
  .note(:upper_limit, "no significant hc signal observed; Bayesian upper limit at 90% C.L. set with systematic uncertainties convolved into likelihood")
  .with_decay_card(decay_card_III)
  .apply(event_sel_III)

algo_III.execute_on([psip_data, psip_incMC, ex_mc_III])