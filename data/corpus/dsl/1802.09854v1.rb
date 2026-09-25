# 1802.09854v1: Study of two-photon decays of pseudoscalar mesons via J/ψ radiative decays
# ψ(3686) → π⁺π⁻ J/ψ, J/ψ → 3γ
# Ordinary analysis (no tagging)

psi3686_data = DatasetManager.real_data.find("709_3686")
psi3686_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for exclusive MC: ψ(3686) → π⁺π⁻ J/ψ, J/ψ → γη → 3γ (representative)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi+ pi- J/psi VPIPI;
    Enddecay

    Decay J/psi
    1.000 gamma eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Psip_to_pipi_Jpsi_eta_3gamma_exclusive_mc"
  config.related_dataset = psi3686_data
  config.events = 500000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

algorithm = Algorithm.new("TwoPhotonPseudoscalar", "00-00-01")

algorithm.set_header(["TwoPhotonPseudoscalarAlg/TwoPhotonPseudoscalar.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .note(:post_fit_cuts, "J/ψ mass window |M(3γ) - m_J/ψ| < 0.08 GeV/c² applied in ROOT")
         .note(:signal_extraction, "π⁰/η/η′/η(1405)/η(1475)/η(1760)/X(1835) signals extracted from M(γγ) using unbinned ML fits in ROOT; sideband background subtraction and angular reweighting of signal MC also performed in ROOT")
         .note(:background_veto, "peaking background from J/ψ → γπ⁰π⁰ estimated via dedicated MC and subtracted in ROOT")
         .with_decay_card(decay_card_signal)

event_selection = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nChrp       "==1"
    nChrn       "==1"
  end
  .assign({chrgp: :pip, chrgn: :pim})
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    15.0
    nGam              ">=3"
  end
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma]) do
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

algorithm.apply(event_selection)
algorithm.execute_on([psi3686_data, psi3686_incMC, exMC_signal])