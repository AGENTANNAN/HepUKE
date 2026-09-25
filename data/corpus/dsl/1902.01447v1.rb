# Search for J/psi -> phi e+ e- at psi(3686)
# Paper: 1902.01447v1, using 448.1e6 psi(3686) events

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 phi e+ e- PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- PHSP;
  Enddecay

  End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_phi_ee_exMC"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("JpsiToPhiEE")
algorithm
  .set_header(["JpsiToPhiEEAlg/JpsiToPhiEE.h"])
  .set_constant({"ECMS" => [:double, 3.686]})

selection = Selection.new
selection.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==3"
  nChrn "==3"
  nNet "==0"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp ">=1"
  nkm ">=1"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :electron, against: [:pion, :kaon]
  nep ">=1"
  nem ">=1"
end
.remove([:kp <= :chrgp, :km <= :chrgn, :ep <= :chrgp, :em <= :chrgn])
.assign({chrgp: :pip, chrgn: :pim})
.kinematic_fit([:pip, :pim, :kp, :km, :ep, :em]) do
  nominal
  constrain_four_momentum
  chi2_cut 40
end

algorithm
  .note(:soft_pion_momentum, "soft pion candidates from psi(3686)->pi+pi-J/psi required to have p < 0.45 GeV/c; no PID applied to pions — identified as remaining tracks after kaon and electron PID removal")
  .note(:ep_cut, "electron candidates required E/p > 0.8 to suppress pion contamination; applied in addition to PID probability requirement")
  .note(:recoil_mass_cut, "recoil mass M(pi+pi-)_rec required in [3.05, 3.15] GeV/c^2 to select J/psi region")
  .note(:best_combination, "if multiple track-to-hypothesis assignments survive PID, the combination with smallest 4C chi2 is retained")
  .note(:signal_region, "signal region defined as |M(pi+pi-)_rec - M_J/psi| < 0.007 GeV/c^2 and |M(K+K-) - M_phi| < 0.010 GeV/c^2; applied in ROOT analysis")
  .note(:sideband_analysis, "sideband method used for background estimation; three types of sideband regions defined around signal region for non-J/psi, non-phi, and non-J/psi-non-phi backgrounds; applied in ROOT analysis")
  .note(:spin_correlation, "J/psi spin correlation from psi(3686) decay and phi helicity angle distribution considered in MC generation; phi->K+K- uses sin^2(theta_K) distribution")
  .note(:photon_conversion_veto, "background from pi0->gamma e+e- conversion and psi(3686)->pi+pi-J/psi with J/psi->phi pi0/eta studied via exclusive MC and found negligible")
  .with_decay_card(decay_card)
  .apply(selection)

algorithm.execute_on([psip_data, psip_incMC, exMC])