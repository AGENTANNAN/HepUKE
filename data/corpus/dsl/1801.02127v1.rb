### Dataset preparation ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

### J/psi -> gamma gamma phi, phi -> K+ K- ###

decay_signal = <<~DECAYCARD
    Decay J/psi
    1.000 gamma gamma phi PHSP;
    Enddecay
    Decay phi
    1.000 K+ K- VSS;
    Enddecay
    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_jpsi_ggphi"
  config.related_dataset = jpsi_data
  config.events = 500000
  config.decay_card = decay_signal
  config.cross_section = :default
end

alg = Algorithm.new("JpsiGammaGammaPhi")
alg.set_header(["JpsiGammaGammaPhiAlg/JpsiGammaGammaPhi.h"])
    .set_constant({"ECMS" => [:double, 3.097]})

sel = Selection.new
sel.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==1"
  nChrn "==1"
  nNet "==0"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp ">=1"
  nkm ">=1"
}
.remove([:kp <= :chrgp, :km <= :chrgn])
.select_photon {
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  tdc_emc_start 0
  tdc_emc_end 14
  nGam ">=2"
}
# Nominal 4C kinematic fit: gamma gamma K+ K-
.kinematic_fit([:gamma, :gamma, :kp, :km]) {
  nominal
  constrain_four_momentum
  chi2_cut 40
}
# Competing hypothesis: 3gamma K+ K- (for veto)
.kinematic_fit([:gamma, :gamma, :gamma, :kp, :km]) {
  constrain_four_momentum
}
# Competing hypothesis: 4gamma K+ K- (for veto)
.kinematic_fit([:gamma, :gamma, :gamma, :gamma, :kp, :km]) {
  constrain_four_momentum
}

alg.note(:background_veto, "veto events where chi2_4c_signal > min(chi2_4c_3gamma, chi2_4c_4gamma) in ROOT")
   .note(:background_veto, "veto pi0: |M(gamma gamma) - m(pi0)| > 0.03 GeV/c^2")
   .note(:background_veto, "veto eta: M(gamma gamma) < 0.50 or > 0.58 GeV/c^2")
   .note(:background_veto, "veto etap: |M(gamma gamma) - m(etap)| > 0.03 GeV/c^2")
   .with_decay_card(decay_signal)
   .apply(sel)

alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])