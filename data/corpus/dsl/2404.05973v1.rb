# Dataset description — energy scan for Ds production
ds_data = DatasetManager.real_data.find("703_4180")
ds_incMC = DatasetManager.inclusive_mc.find("703_4180")

# Decay card: Ds+ → pi+ phi, phi → e+ e- (via D_s*+ D_s- production)
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s- PHSP;
  Enddecay
  Decay D_s*+
  1.0000 D_s+ gamma VSP_PWAVE;
  Enddecay
  Decay D_s+
  1.0000 pi+ phi SVS;
  Enddecay
  Decay D_s-
  1.0000 K+ K- pi- PHSP;
  Enddecay
  Decay phi
  1.0000 e+ e- VLL;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dsp_to_pip_phi_ee"
  config.related_dataset = ds_data
  config.events = 500000
  config.decay_card = decay_card
  config.cross_section = :default
end

# Algorithm — Mode I: Ds+ → pi+ phi, phi → e+ e-
alg_mode1 = Algorithm.new("DspToPiPhiEE")
alg_mode1.set_header(["DspToPiPhiEEAlg/DspToPiPhiEE.h"])
alg_mode1.set_constant({"ECMS" => [:double, 4.178]})

sel_mode1 = Selection.new
sel_mode1.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=1"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=1"
}
.pid(method: :probability) {
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  identify :pion, against: [:kaon, :proton]
  identify :kaon, against: [:pion, :proton]
}

alg_mode1
  .note(:electron_pid_criteria, "Electron PID: L'(e) > 0.001 and L'(e)/(L'(e)+L'(pi)+L'(K)) > 0.8; E/p > 0.8 (0.7) for p > (<) 0.4 GeV/c")
  .note(:gamma_conversion_veto, "Veto e+e- pairs from gamma conversions: discard events with e+e- vertex distance from IP in [2.0, 8.0] cm; photon recovery within 5-degree cone around e+-")
  .note(:ds_mass_window, "Ds+ candidate invariant mass required in [1.88, 2.02] GeV/c^2")
  .note(:recoil_mass_deltaM, "Signal candidates selected from M_rec vs DeltaM 2D plane; M_rec and DeltaM windows optimized per energy point via S/sqrt(S+B)")
  .note(:phi_mass_window, "M(e+e-) required in [0.98, 1.04] GeV/c^2 (phi mass window)")
  .note(:ds_star_photon, "D_s*+ → gamma D_s+ transition photon: if multiple photons, choose one with recoil mass of D_s+ gamma closest to m_Ds*+")
  .with_decay_card(decay_card)
  .apply(sel_mode1)

root_files = alg_mode1.execute_on([ds_data, ds_incMC, exMC_signal])