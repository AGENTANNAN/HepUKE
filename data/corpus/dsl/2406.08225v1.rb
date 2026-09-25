# Paper: 2406.08225v1
# ψ(3686) → γ 2(π+π-)η via radiative transitions
# Observation of ηc(1S,2S) and χcJ decays to 2(π+π-)η

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: ψ(3686) → γ ηc, ηc → 2(π+π-)η, η → γγ
# All intermediate states (ηc(1S), χcJ, ηc(2S)) share the same final state
# and selection criteria → single Algorithm
decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 pi+ pi+ pi- pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_gamma_2pip2pim_eta"
  config.related_dataset = psip_data
  config.events = 500000
  config.decay_card = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("PsipGamma2PiPiEta")
algorithm.set_header(["PsipGamma2PiPiEtaAlg/PsipGamma2PiPiEta.h"])
          .set_constant({ "ECMS" => [:double, 3.686] })

# Event selection: 4 charged tracks (all pions), ≥3 photons
# η → γγ reconstructed via Kalman fit, then 4C kinematic fit
# with the radiative γ + 4π + η
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
  nGam ">=3"
end
.assign({ chrgp: :pip, chrgn: :pim })
# Reconstruct η → γγ with 1C Kalman fit
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 20
  neta ">=1"
end
# 4C kinematic fit: radiative γ + 4π + η under 4-momentum conservation
.kinematic_fit([:gamma, :pip, :pip, :pim, :pim, :eta]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

algorithm
  .note(:competing_hypothesis_veto,
    "Events with χ²_4C(3γ 2(π+π-)) > χ²_4C(4γ 2(π+π-)) vetoed to suppress background with extra photons;
     Events with χ²_4C(3γ 2(π+π-)) > χ²_4C(2γ 2(π+π-)) vetoed to suppress background with missing photons")
  .note(:chi2_5c_cut,
    "Optimal χ²_5C < 15 applied in ROOT after s/√(s+b) optimisation")
  .note(:background_veto,
    "ψ(3686)→P J/ψ vetoed by M_recoil(P) outside J/ψ mass window;
     ψ(3686)→π0H vetoed by M(γγ) outside π0 mass region;
     For γηc(2S)/χcJ modes, ψ(3686)→γχcJ, χcJ→γJ/ψ vetoed")
  .with_decay_card(decay_card)
  .apply(event_selection)

algorithm.execute_on([psip_data, psip_incMC, exMC_signal])