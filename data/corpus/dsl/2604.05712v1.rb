# CKM angle gamma with a novel unbinned optimal Fourier method (Letter version).
# BESIII: quantum-correlated D Dbar decays at psi(3770), 8 fb^-1 (2010-2011 + 2021-2022).
# Signal D -> Ks/L h'+ h'- decays with three tag categories: flavor / CP-eigenstate / self-conjugate.
# LHCb part is out of BOSS scope.

### Datasets ###
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Signal decay cards: D0 -> Ks pi+ pi- and D0 -> Ks K+ K- (D_DALITZ dynamical model).
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

### Event selection (BOSS) — two DT algorithms (Rule T1) ###

# --- Channel A: D0 -> Ks pi+ pi- signal, DT with all tag categories ---
alg_A_name = "DKsPiPiDTagPRL"
alg_A = TagAnalysis.new(alg_A_name)
alg_A.set_header(["#{alg_A_name}Alg/#{alg_A_name}.h"])
     .set_constant({ "ECMS" => [:double, 3.773] })
     .set_alias({ "std::vector<double>" => "Vdouble" })

alg_A.tag_side(:D0) do |t|
  t.modes :D0toKsPiPi
  t.charm 1
end

alg_A.tag_side(:D0) do |t|
  # flavor-specific + CP-even + CP-odd + self-conjugate Ks h h tags
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi,
          :D0toKK, :D0toPiPi, :D0toKsPi0Pi0, :D0toPiPiPi0,
          :D0toKsPi0, :D0toKsEta, :D0toKsEtaPr, :D0toKsOmega,
          :D0toKsPiPi, :D0toKsKK
  t.charm -1
  t.rank_by :inv
end

alg_A.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_A.note(:mass_constraints,
           "D, K_S0, K_L0 candidates are constrained to their nominal masses to mitigate Dalitz-plot resolution effects.")
     .note(:selection_source,
           "Final-state particle selection follows the previous BESIII binned strong-phase measurement (JHEP 06 (2025) 086).")
     .note(:partial_reconstruction,
           "For K_L0 h+ h- signal (and K_L0 in tag), the K_L0 is treated as a missing particle: M_miss^2 = E_miss^2 - |p_miss|^2 is used as the fit variable at the ROOT stage.")
     .note(:optimal_fourier_weights,
           "Per-event optimal Fourier weights (Ref. [26] method) with M_pi=2 and M_K=1 are applied at the ROOT stage in the joint BESIII+LHCb fit for gamma.")
     .note(:peaking_backgrounds,
           "Dominant peaking backgrounds are D -> pi+ pi- h+ h- for Ks h+h- signal and D -> Ks(->pi0 pi0) h+ h- for KL h+ h- signal; handled in the ROOT fit.")

alg_A.apply
alg_A.execute_on([psi3770_data, psi3770_incMC, exMC_KsPiPi])

# --- Channel B: D0 -> Ks K+ K- signal, DT with all tag categories ---
alg_B_name = "DKsKKDTagPRL"
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
          :D0toKsPi0, :D0toKsEta, :D0toKsEtaPr, :D0toKsOmega,
          :D0toKsPiPi, :D0toKsKK
  t.charm -1
  t.rank_by :inv
end

alg_B.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_B.note(:mass_constraints,
           "D, K_S0, K_L0 candidates are constrained to their nominal masses.")
     .note(:peaking_backgrounds_kk,
           "Dominant peaking backgrounds D -> pi+pi-K+K- (~0.5%) and D -> Ks(->pi0 pi0)K+K- (~2.5%); subtracted in ROOT.")
     .note(:optimal_fourier_weights,
           "Per-event optimal Fourier weights (M_pi=2, M_K=1) applied at ROOT stage; joint fit with LHCb determines gamma.")

alg_B.apply
alg_B.execute_on([psi3770_data, psi3770_incMC, exMC_KsKK])
