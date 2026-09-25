# Paper 2108.02405v3: D0 -> omega phi at psi(3770)
# omega -> pi+ pi- pi0, phi -> K+ K-, pi0 -> gamma gamma
# Single-tag technique at sqrt(s)=3.773 GeV

decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 VSS;
  Enddecay
  Decay D0
  1.0000 omega phi PHSP;
  Enddecay
  Decay omega
  1.0000 pi+ pi- pi0 PHSP;
  Enddecay
  Decay phi
  1.0000 K+ K- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

sigMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "sig_D0_omega_phi"
  config.related_dataset = data_3773
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("D0OmegaPhi")
algorithm.set_header(["D0OmegaPhiAlg/D0OmegaPhi.h"])
algorithm.set_constant({"ECMS" => [:double, 3.773]})

event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=2"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :pion, against: [:kaon]
  identify :kaon, against: [:pion]
  nkp ">=1"
  nkm ">=1"
  npip ">=1"
  npim ">=1"
end
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
end
.kinematic_fit([:pip, :pim, :pi0, :kp, :km]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

algorithm
  .note(:single_tag_technique, "single-tag technique: only one D0 is fully reconstructed via omega phi; D0_bar decays generically")
  .note(:pi0_reconstruction, "pi0 -> gamma gamma with Kalman kinematic fit constraining M(gamma gamma) to nominal pi0 mass; mass window (0.115,0.150) GeV/c^2 before fit")
  .note(:pid_method, "PID uses likelihood comparison: L_pi > L_K for pi+-, L_K > L_pi for K+-")
  .note(:d0_selection, "D0 signal: M_BC > 1.84 GeV/c^2; -0.03 < DeltaE < 0.02 GeV; best candidate: minimum |DeltaE|")
  .note(:ks_veto, "K_S0 veto: remove events with M(pi+pi-) in (0.490,0.503) GeV/c^2")
  .note(:signal_extraction, "2D unbinned maximum likelihood fit to M(pi+pi-pi0) vs M(K+K-) in M_BC signal and sideband regions; SIGNAL, BKGI, BKGII, BKGIII components; sideband subtraction for non-D0 backgrounds")
  .note(:polarization, "angular analysis of cos(theta_omega) and cos(theta_K) distributions for longitudinal polarization fraction f_L measurement; efficiency-corrected yields in 5 equal bins")
  .with_decay_card(decay_card)
  .apply(event_selection)

algorithm.execute_on([data_3773, incMC_3773, sigMC])