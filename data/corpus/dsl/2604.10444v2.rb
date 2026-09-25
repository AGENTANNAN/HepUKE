# Amplitude analysis and BF measurement of D+ -> pi+ pi+ pi- eta and D+ -> pi+ pi0 pi0 eta.
# BESIII, psi(3770), 20.3 fb^-1. Double-tag (DT) technique.
# Signal side: D+ -> pi+ pi+(0) pi-(0) eta ; Tag side: D- reconstructed in 5 hadronic modes.
# Rule T1: two independent signal channels => two TagAnalysis algorithms.

### Datasets ###
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for D+ -> pi+ pi+ pi- eta (three-pion + eta) signal; D- -> K+ pi- pi- tag
decay_card_3pi_eta = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D+   D-                            VSS;
    Enddecay

    Decay D+
    1.0000  pi+  pi+  pi-  eta                 PHSP;
    Enddecay

    Decay D-
    1.0000  K+   pi-  pi-                      PHSP;
    Enddecay

    Decay eta
    1.0000  gamma  gamma                       PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for D+ -> pi+ pi0 pi0 eta signal
decay_card_pi_pi0pi0_eta = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D+   D-                            VSS;
    Enddecay

    Decay D+
    1.0000  pi+  pi0  pi0  eta                 PHSP;
    Enddecay

    Decay D-
    1.0000  K+   pi-  pi-                      PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                       PHSP;
    Enddecay

    Decay eta
    1.0000  gamma  gamma                       PHSP;
    Enddecay

    End
DECAYCARD

exMC_3pi_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dp_3pi_eta_signal"
  config.related_dataset = psi3770_data
  config.events          = 500000
  config.decay_card      = decay_card_3pi_eta
  config.cross_section   = :default
end
exMC_3pi_eta.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_pi_pi0pi0_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dp_pi_pi0pi0_eta_signal"
  config.related_dataset = psi3770_data
  config.events          = 500000
  config.decay_card      = decay_card_pi_pi0pi0_eta
  config.cross_section   = :default
end
exMC_pi_pi0pi0_eta.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) — Double-tag D+ analysis ###

# --- Signal channel A: D+ -> pi+ pi+ pi- eta with eta -> gamma gamma ---
alg_A_name = "Dp_PiPiPiEta"
alg_A = TagAnalysis.new(alg_A_name)
alg_A.set_header(["#{alg_A_name}Alg/#{alg_A_name}.h"])
     .set_constant({ "ECMS" => [:double, 3.773] })
     .set_alias({ "std::vector<double>" => "Vdouble" })

# Tag side: D- in 5 hadronic modes
alg_A.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

# Signal side: 3 charged pions (net charge +1) + 2 photons for eta
alg_A.signal_side do |s|
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_A.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 100
end

alg_A.note(:signal_region,
           "Signal region 1.863 < MBC < 1.877 GeV/c^2; DeltaE(sig) in (-0.033, 0.032) GeV; DeltaE(tag) mode-dependent windows applied in ROOT.")
     .note(:ks0_veto,
           "Ks0 veto: M(pi+ pi-) not in (0.478, 0.517) GeV/c^2 to reject D+ -> Ks0(->pi+pi-) pi+ eta.")
     .note(:etap_veto,
           "eta' veto: M(pi+ pi- eta) not in (0.8, 1.0) GeV/c^2 to reject D+ -> eta'(->pi+pi-eta) pi+.")
     .note(:best_candidate,
           "Best signal candidate chosen by minimum |DeltaE_sig|; best tag by minimum |DeltaE_tag|.")
     .note(:amplitude_analysis,
           "Amplitude analysis performed with an additional D+ mass constraint in a second kinematic fit; unbinned maximum-likelihood method with isobar model (BW / Gounaris-Sakurai / Flatte) at ROOT stage.")

alg_A.apply
alg_A.execute_on([psi3770_data, psi3770_incMC, exMC_3pi_eta])

# --- Signal channel B: D+ -> pi+ pi0 pi0 eta with eta -> gamma gamma and pi0 -> gamma gamma ---
alg_B_name = "Dp_PiPi0Pi0Eta"
alg_B = TagAnalysis.new(alg_B_name)
alg_B.set_header(["#{alg_B_name}Alg/#{alg_B_name}.h"])
     .set_constant({ "ECMS" => [:double, 3.773] })
     .set_alias({ "std::vector<double>" => "Vdouble" })

alg_B.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

# Signal side: 1 charged pion (charge +1) + 6 photons (2 pi0 + eta -> gamma gamma)
alg_B.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 1
  s.photons 6
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_B.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 50
end

alg_B.note(:signal_region,
           "Signal region 1.863 < MBC < 1.877 GeV/c^2; DeltaE(sig) in (-0.052, 0.038) GeV; DeltaE(tag) mode-dependent windows applied in ROOT.")
     .note(:ks0_veto_pi0pi0,
           "Ks0 veto on pi0 pi0 pair: M(pi0 pi0) not in (0.448, 0.534) GeV/c^2.")
     .note(:etap_veto,
           "eta' veto on pi0(0) eta invariant masses not in (0.8, 1.0) GeV/c^2.")
     .note(:pi_4pi0_veto,
           "Reject events containing an accompanying D+ -> pi+ pi0 pi0 pi0 combination with DeltaE in (-0.1, 0.1) GeV and MBC in (1.83, 1.89) GeV/c^2 to suppress this peaking background (~77% background rejected, ~98% signal retained).")
     .note(:best_candidate,
           "Best candidate chosen by minimum |DeltaE_sig|.")
     .note(:amplitude_analysis,
           "Amplitude analysis performed with an additional D+ mass constraint; isobar-model unbinned maximum-likelihood at ROOT stage.")

alg_B.apply
alg_B.execute_on([psi3770_data, psi3770_incMC, exMC_pi_pi0pi0_eta])
