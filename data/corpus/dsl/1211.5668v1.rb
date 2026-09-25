# =============================================================================
# BESIII: J/psi -> gamma omega phi (X(1810) PWA); omega -> pi+ pi- pi0, phi -> K+K-
# arXiv:1211.5668v1
# =============================================================================

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card_gammaX1810 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma X(1810) PHSP;
  Enddecay

  Decay X(1810)
  1.0000 omega phi PHSP;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay phi
  1.0000 K+ K- VSS;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_gwphi = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_gamma_omega_phi"
  c.related_dataset = jpsi_data
  c.events          = 500000
  c.decay_card      = decay_card_gammaX1810
  c.cross_section   = :default
end

alg = Algorithm.new("JpsiGammaOmegaPhi")
alg.set_header(["JpsiGammaOmegaPhiAlg/JpsiGammaOmegaPhi.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      Vz  10.0    # |Vz| < 10 cm
      Vr  1.0     # |Vr| < 1 cm
      nChrp "==2"
      nChrn "==2"
      nNet  "==0"
    }
    .select_photon {
      nGam ">=3"
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      angle_to_track    10.0
      tdc_emc_start     0
      tdc_emc_end       14
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      identify :pion, against: [:kaon, :proton]
      nkp  "==1"
      nkm  "==1"
      npip "==1"
      npim "==1"
    }
    # Nominal 4C fit: J/psi -> 3 gamma K+ K- pi+ pi-
    .kinematic_fit([:gamma, :gamma, :gamma, :kp, :km, :pip, :pim]) {
      nominal
      constrain_four_momentum
      chi2_cut 200   # loose; tight cut (40) applied in ROOT after S/sqrt(S+B) optim.
    }
    # Competing-hypothesis fit: 2 gamma K+K- pi+pi- (no chi2_cut, no nominal)
    .kinematic_fit([:gamma, :gamma, :kp, :km, :pip, :pim]) {
      constrain_four_momentum
    }
    # Competing-hypothesis fit: 4 gamma K+K- pi+pi- (no chi2_cut, no nominal)
    .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :kp, :km, :pip, :pim]) {
      constrain_four_momentum
    }

alg.note(:pi0_mass_window,
         "reconstruct pi0 from the gamma-gamma pair (of the three selected photons) "\
         "whose invariant mass is closest to the nominal pi0 mass; require "\
         "|M(gamma gamma) - M_pi0| < 20 MeV/c^2")
   .note(:omega_signal_window,
         "|M(pi+ pi- pi0) - M_omega| < 40 MeV/c^2 for the omega signal region")
   .note(:phi_signal_window,
         "|M(K+ K-) - M_phi| < 15 MeV/c^2 for the phi signal region")
   .note(:etap_background_veto,
         "require M(gamma pi+ pi- pi0) > 1.0 GeV/c^2 to remove J/psi -> phi eta' "\
         "(eta' -> gamma omega) background")

alg.with_decay_card(decay_card_gammaX1810).apply(sel)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_gwphi])
