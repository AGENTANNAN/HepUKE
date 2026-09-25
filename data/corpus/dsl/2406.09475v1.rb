# Paper: 2406.09475v1
# J/ψ → ω K+K-η, search for X(1870) → K+K-η
# ω → π+π-π0, π0 → γγ, η → γγ

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for J/ψ → ω K+K-η (phase space)
decay_card_omega_kpi_eta = <<~DECAYCARD
    Decay J/psi
    1.0000 omega K+ K- eta PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for J/ψ → ω X(1870), X(1870)→K+K-η
decay_card_x1870 = <<~DECAYCARD
    Decay J/psi
    1.0000 omega X_1870 PHSP;
    Enddecay

    Decay X_1870
    1.0000 K+ K- eta PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_omega_kpi_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_omega_kpi_eta"
  config.related_dataset = jpsi_data
  config.events = 500000
  config.decay_card = decay_card_omega_kpi_eta
  config.cross_section = :default
end

exMC_x1870 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_omega_x1870_kpi_eta"
  config.related_dataset = jpsi_data
  config.events = 500000
  config.decay_card = decay_card_x1870
  config.cross_section = :default
end

algorithm = Algorithm.new("JpsiOmegaKKEta")
algorithm.set_header(["JpsiOmegaKKEtaAlg/JpsiOmegaKKEta.h"])
          .set_constant({ "ECMS" => [:double, 3.097] })

# Event selection: K+K-π+π-4γ, zero net charge
# ω → π+π-π0(→γγ), η → γγ
# 5C kinematic fit: 4-momentum conservation + one γγ to π0 mass
event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=4"
end
.assign({ chrgp: :pip, chrgn: :pim })
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  identify :pion, against: [:kaon, :proton]
end
# 5C kinematic fit: K+ K- π+ π- + 4γ with 4-momentum + one γγ→π0 constraint
.kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma, :gamma]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 80
end

algorithm
  .note(:d0_vz_cuts,
    "Rxy < 1.0 cm, |Vz| < 10.0 cm; cosθ(MDC) < 0.93")
  .note(:mass_windows,
    "|M(π+π-π0) - m(ω)| < 0.02 GeV/c² and |M(γγ) - m(η)| < 0.02 GeV/c² applied in ROOT after kinematic fit")
  .note(:etap_veto,
    "|M(π+π-η) - m(η')| > 0.025 GeV/c² applied in ROOT to suppress η' background")
  .note(:x1870_search,
    "Search for X(1870) → K+K-η via M(K+K-η) distribution; BW signal shape with mass-dependent width")
  .note(:omega_dalitz,
    "ω simulated with Dalitz plot distribution; X(1870) simulated with Breit-Wigner")
  .with_decay_card(decay_card_omega_kpi_eta)
  .apply(event_selection)

algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_omega_kpi_eta, exMC_x1870])