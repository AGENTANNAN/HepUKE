# ============================================================
# Dataset preparation
# ============================================================
# psi(3770) data (20.3 fb^-1 at sqrt(s) = 3.773 GeV)
psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# --------------------------------------------------------------------
# Signal decay card: D+ -> K_S0 K_L0 pi+
# --------------------------------------------------------------------
decay_card_signal = <<~DECAYCARD
  Decay D+
  1.0000  K_S0  K_L0  pi+   PHSP;
  Enddecay
  CDecay D-

  Decay K_S0
  1.0000  pi+  pi-           PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Dp_to_KsKLPi"
  config.related_dataset = psi3770_data
  config.events         = 500000
  config.decay_card     = decay_card_signal
  config.cross_section  = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

# ============================================================
# Event selection: DTag double tag at psi(3770).
# ST: D- reconstructed in 3 hadronic tag modes (K+ pi- pi-, K_S0 pi-, K_S0 pi+ pi- pi-)
# DT: D+ -> K_S0 K_L0 pi+ from unused tracks; K_L0 treated as missing particle
# ============================================================

alg = TagAnalysis.new("DpToKsKLPi")
alg.set_header(["DpToKsKLPiAlg/DpToKsKLPi.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })
   .with_decay_card(decay_card_signal)

# Single-tag D- reconstructed in three hadronic modes
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKsPiPiPi
  t.charm -1
end

# Signal side: K_S0(->pi+ pi-) + pi+ (3 charged tracks total: 2 pi+, 1 pi-)
# K_L0 is missing; net charge = +1
alg.signal_side do |s|
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
  s.missing :K_L0
  s.min_photon_angle 10.0
end

# Kinematic fit: total 4-momentum, D+/- mass constraint, K_S0 mass constraint;
# missing 4-momentum interpreted as K_L0; require chi^2 < 100.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 100
end

alg.note(:ks_reconstruction,
         "K_S0 reconstructed via K_S0 -> pi+ pi- with secondary vertex fit, using the standard BESIII selection " \
         "(Refs. [35,36]); at least one oppositely-charged pion pair consistent with K_S0.")
   .note(:pi0_eta_veto,
         "Veto events with a reconstructed pi0 or eta -> gamma gamma to suppress " \
         "D+ -> K_S0 pi+ eta (peaking in Mmiss^2 at m_eta) and D+ -> K_S0 K_S0 pi+ with one K_S0 -> pi0 pi0.")
   .note(:kl_missing_mass,
         "Missing mass squared Mmiss^2 required in (0.220, 0.265) GeV^2/c^4 (K_L0 peak). " \
         "Applied in ROOT stage after kinematic fit.")
   .note(:tag_side,
         "ST D- reconstructed in K+ pi- pi-, K_S0 pi-, K_S0 pi+ pi- pi-; " \
         "M_BC in (1.863, 1.877) GeV/c^2; best candidate per event chosen by min |Delta E|.")

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_signal])
