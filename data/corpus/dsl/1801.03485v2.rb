### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

### Decay cards ###
# Signal: psi(2S) -> pi0 hc, hc -> pi+ pi- J/psi, J/psi -> l+ l-
decay_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 hc PHSP;
    Enddecay
    Decay hc
    1.000 pi+ pi- J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- VLL;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Reference: psi(2S) -> eta J/psi, eta -> pi0 pi+ pi-, J/psi -> l+ l-
decay_ref = <<~DECAYCARD
    Decay psi(2S)
    1.000 eta J/psi PHSP;
    Enddecay
    Decay eta
    1.000 pi0 pi+ pi- PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- VLL;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_hc_pipi_jpsi"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_signal
  config.cross_section = :default
end

exMC_ref = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "ref_eta_jpsi_pi0pipi"
  config.related_dataset = psip_data
  config.events = 500000
  config.decay_card = decay_ref
  config.cross_section = :default
end

### Event selection ###
alg = Algorithm.new("hc_pipi_jpsi")
alg.set_header(["hc_pipi_jpsiAlg/hc_pipi_jpsi.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

sel = Selection.new
sel.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
  nNet "==0"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0, treat_as_electron_if_energy_above: 0.6
  nlp ">=1"
  nlm ">=1"
}
.remove([:lp <= :chrgp, :lm <= :chrgn])
.assign({chrgp: :pip, chrgn: :pim})
.select_photon {
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  tdc_emc_start 0
  tdc_emc_end 14
  nGam ">=2"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=1"
}
# 5C kinematic fit: 4C + pi0 mass constraint (from prior kalman fit)
.kinematic_fit([:pi0, :pip, :pim, :lp, :lm]) {
  nominal
  constrain_four_momentum
  chi2_cut 60
}

alg.note(:pid_correction_method, "lepton ID: E/pc > 0.7 -> electron, E/pc < 0.3 -> muon; pions assigned by p < 1 GeV/c")
   .note(:background_veto, "M(pi+pi-) > 0.3 GeV/c^2 to reject gamma conversion background")
   .note(:efficiency_curve, "pi0 mass window: 120 < M(gamma gamma) < 145 MeV/c^2")
   .with_decay_card(decay_signal)
   .apply(sel)

alg.execute_on([psip_data, psip_incMC, exMC_signal, exMC_ref])