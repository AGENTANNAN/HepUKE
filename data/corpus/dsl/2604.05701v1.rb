# BESIII portion of the joint BESIII+LHCb CKM angle gamma measurement.
# Quantum-correlated D Dbar decays at psi(3770), 8 fb^-1 (2010-2011 + 2021-2022).
# Signal D -> Ks/L h'+ h'- decays, tagged against flavor, CP-even, CP-odd, and other Ks h h tags.

### Dataset ###
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for signal D -> Ks pi+ pi- with anti-D decaying to CP-even/flavor tag(s)
decay_card_KsPiPi = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D0    anti-D0                     VSS_MIX;
    Enddecay

    Decay D0
    1.0000  K_S0  pi+  pi-                    D_DALITZ;
    Enddecay

    Decay anti-D0
    1.0000  K+  pi-                           PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                          PHSP;
    Enddecay

    End
DECAYCARD

decay_card_KsKK = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D0    anti-D0                     VSS_MIX;
    Enddecay

    Decay D0
    1.0000  K_S0  K+  K-                      D_DALITZ;
    Enddecay

    Decay anti-D0
    1.0000  K+  pi-                           PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                          PHSP;
    Enddecay

    End
DECAYCARD

exMC_KsPiPi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D_KsPiPi_DT_signal"
  config.related_dataset = psi3770_data
  config.events          = 500000
  config.decay_card      = decay_card_KsPiPi
  config.cross_section   = :default
end
exMC_KsPiPi.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_KsKK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D_KsKK_DT_signal"
  config.related_dataset = psi3770_data
  config.events          = 500000
  config.decay_card      = decay_card_KsKK
  config.cross_section   = :default
end
exMC_KsKK.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
# Rule T1: two independent signal decay modes => two TagAnalysis algorithms.

# --- Signal channel A: D -> Ks pi+ pi- vs flavor / CP tag ---
alg_A_name = "DKsPiPiDTag"
alg_A = TagAnalysis.new(alg_A_name)
alg_A.set_header(["#{alg_A_name}Alg/#{alg_A_name}.h"])
     .set_constant({ "ECMS" => [:double, 3.773] })
     .set_alias({ "std::vector<double>" => "Vdouble" })

# Double-tag: signal D0 -> Ks pi+ pi- (via :D0toKsPiPi) paired with the many tag modes.
alg_A.tag_side(:D0) do |t|
  t.modes :D0toKsPiPi
  t.charm 1
end

# Tag side: 11 hadronic + CP + Ks h h tag modes covering flavor / CP-even / CP-odd
alg_A.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi,          # flavor tags K pi (pi0/pi+pi-)
          :D0toKK, :D0toPiPi, :D0toKsPi0Pi0, :D0toPiPiPi0, # CP-even tags
          :D0toKsPi0, :D0toKsEta, :D0toKsEtaPr,          # CP-odd tags
          :D0toKsOmega, :D0toKsPiPi, :D0toKsKK            # other Ks h h tags (self-conjugate DT)
  t.charm -1
  t.rank_by :inv
end

alg_A.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_A.note(:signal_channel, "Signal D0 -> Ks0 pi+ pi- decay, reconstructed via DTag Ks pi+ pi- with the tag side chosen from three categories: flavor (K pi, K pi pi0, K pi pi pi, K e nu_e), CP-even (KK, pi pi, Ks pi0 pi0, pi pi pi0), CP-odd (Ks pi0, Ks eta(->gg or pi pi pi0), Ks eta'(->pi pi gamma or pi pi eta), Ks omega(->pi pi pi0)), and other self-conjugate Ks h+ h- modes.")
   .note(:mass_constraints,
         "D, K_S0, K_L0 candidates are constrained to their PDG masses to mitigate Dalitz-plot resolution effects.")
   .note(:missing_kl,
         "For partial reconstruction K_L0 pi+ pi- and K_L0 K+ K- tags, missing mass squared M_miss^2 is used as the fit variable; for D -> K e nu_e the observable is U_miss = E_miss - |p_miss|.")
   .note(:signal_yield_extraction,
         "1D unbinned maximum-likelihood fits on M_BC of each mode extract signal yields; double-sided Crystal Ball convolved with Gaussian; background from inclusive MC.")
   .note(:dalitz_efficiency,
         "Dalitz-plot efficiency map fitted by 2D polynomial (up to 4th order) from uniform PHSP signal MC.")
   .note(:optimal_fourier_weights,
         "Per-event optimal Fourier weights (Ref. [21] method) applied at the ROOT stage.")

alg_A.apply
alg_A.execute_on([psi3770_data, psi3770_incMC, exMC_KsPiPi])

# --- Signal channel B: D -> Ks K+ K- vs flavor / CP tag ---
alg_B_name = "DKsKKDTag"
alg_B = TagAnalysis.new(alg_B_name)
alg_B.set_header(["#{alg_B_name}Alg/#{alg_B_name}.h"])
     .set_constant({ "ECMS" => [:double, 3.773] })
     .set_alias({ "std::vector<double>" => "Vdouble" })

alg_B.tag_side(:D0) do |t|
  t.modes :D0toKsKK
  t.charm 1
end

alg_B.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi,
          :D0toKK, :D0toPiPi, :D0toKsPi0Pi0, :D0toPiPiPi0,
          :D0toKsPi0, :D0toKsEta, :D0toKsEtaPr,
          :D0toKsOmega, :D0toKsPiPi, :D0toKsKK
  t.charm -1
  t.rank_by :inv
end

alg_B.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_B.note(:signal_channel, "Signal D0 -> Ks0 K+ K- decay, reconstructed via DTag Ks K K with the tag side chosen from the same flavor / CP-even / CP-odd / Ks h h categories as the Ks pi pi channel.")
     .note(:peaking_background,
           "Peaking backgrounds D -> pi+pi-K+K- (0.5%) and D -> Ks(->pi0 pi0)K+K- (2.5%) are handled at the ROOT fit stage.")
     .note(:mass_constraints,
           "D, K_S0, K_L0 candidates constrained to PDG masses.")
     .note(:optimal_fourier_weights,
           "Per-event optimal Fourier weights (Ref. [21]) applied at the ROOT stage in the joint BESIII+LHCb fit for the CKM angle gamma.")

alg_B.apply
alg_B.execute_on([psi3770_data, psi3770_incMC, exMC_KsKK])
