# Improved measurement of absolute BF of inclusive Λc+ → K_S^0 X
# Paper: 2502.20821v2
# DT method at 7 energy points (4599.53–4698.82 MeV).
# Tag side: anti-Λc- reconstructed via 11 hadronic tag modes (8 available).
# Signal side: K_S^0 → π+π- from remaining tracks recoiling against ST anti-Λc-.

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

# Seven centre-of-mass energy points (Table 2)
data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

scan_data  = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]
scan_incMC = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

decay_card_signal = <<~DECAYCARD
  Decay vpho
  1 Lambda_c+ anti-Lambda_c- J2BB1;
  Enddecay
  Decay Lambda_c+
  1.0000 K_S0 X PHSP;
  Enddecay
  Decay anti-Lambda_c-
  1.0000 anti-p- K_S0    PHSP;
  1.0000 anti-p- K+ pi-  PHSP;
  1.0000 anti-p- K_S0 pi0  PHSP;
  1.0000 anti-p- K_S0 pi- pi+  PHSP;
  1.0000 anti-p- K+ pi- pi0  PHSP;
  1.0000 anti-Lambda pi-  PHSP;
  1.0000 anti-Lambda pi- pi0  PHSP;
  1.0000 anti-Lambda pi- pi+ pi-  PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  Decay anti-Lambda
  1.0000 anti-p- pi+ PHSP;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "sig_lambdac_ksx"
  config.events        = 1_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

analysis = TagAnalysis.new("LambdacKsXInclusive")
analysis
  .set_header(["LambdacKsXInclusiveAlg/LambdacKsXInclusive.h"])

# ST anti-Λc- tag: 8 of 11 modes have authoritative symbols.
# The 3 modes with Σ0/Σ- in the final state have no matching DTagAlg symbol
# and are dropped.
analysis.tag_side(:Lambdac) do
  modes :LambdacPtoKsP, :LambdacPtoKPiP, :LambdacPtoKsPi0P,
        :LambdacPtoKsPiPiP, :LambdacPtoKPiPi0P,
        :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
        :LambdacPtoLambdaPiPiPi
  charm -1
end

# Signal side: K_S^0 → π+π- from tracks not used by the tag
analysis.signal_side do
  charged [:piplus, :piminus]
end

analysis.fit do
  constrain_four_momentum
  chi2_cut 200
end

analysis
  .note(:tag_modes_dropped,
    "3 of 11 tag modes from the paper involve Σ0/Σ- (anti-Σ0 π-, anti-Σ- π0, anti-Σ- π-π+) " \
    "and have no matching DTagAlg symbol in the authoritative vocabulary. These are dropped. " \
    "The 8 retained modes are: anti-p K_S^0, anti-p K+π-, anti-p K_S^0 π0, anti-p K_S^0 π-π+, " \
    "anti-p K+π-π0, anti-Λ π-, anti-Λ π-π0, anti-Λ π-π+π-.")
  .note(:ks_reconstruction,
    "K_S^0 reconstructed from two oppositely charged tracks with: " \
    "|M(π+π-)| ∈ [487, 511] MeV/c^2 (~3σ mass window), " \
    "vertex fit χ² < 100, decay length L > 2σ_L. " \
    "If multiple K_S^0 candidates, the one with max L/σ_L is selected.")
  .note(:antiLambda_veto,
    "For tag modes anti-p K_S^0 π0, anti-p K_S^0 π-π+, and anti-Σ- π-π+: " \
    "anti-Λ veto via M(anti-p π+) ∉ (1110, 1120) MeV/c^2. " \
    "For anti-p K_S^0 π0: anti-Σ- veto via M(anti-p π0) ∉ (1170, 1200) MeV/c^2. " \
    "For anti-Λ π-π+π-, anti-Σ- π0, anti-Σ- π-π+: K_S^0 veto via M(π+π-), M(π0π0) ∉ (480, 520) MeV/c^2.")
  .note(:st_selection,
    "ST anti-Λc- yields extracted via unbinned ML fit to M_BC distribution. " \
    "Signal shape: MC-simulated convolved with Gaussian (data-MC resolution difference). " \
    "Background shape: ARGUS function. ΔE windows are mode-dependent (Table 3). " \
    "If multiple candidates per mode, the one with min |ΔE| is chosen.")
  .note(:dt_fit,
    "DT signal yields extracted via simultaneous 2D unbinned ML fit to M_BC vs M(π+π-) " \
    "across all 7 energy points. Signal shapes from matched MC events convolved with Gaussians. " \
    "Background modeled with Chebyshev polynomials. This fit is implemented entirely in ROOT.")
  .note(:ks_efficiency_correction,
    "K_S^0 reconstruction efficiency correction factor (1.8 ± 0.3)% determined from " \
    "control sample J/ψ → K_S^0 K± π∓, applied as a function of K_S^0 momentum.")
  .note(:rdp_match,
    "Truth-matching via R_dp = |p_truth − p_rec| / |p_truth| < 0.5 for K_S^0 daughter π± " \
    "used to separate matched vs unmatched events for signal shape determination.")
  .note(:bf_formula,
    "B(Λc+ → K_S^0 X) treated as shared parameter across energy points in simultaneous fit: " \
    "N_DT_α → B(K_S^0→π+π-) · Σ_i (N_ST_iα / ε_ST_iα · ε_DT_iα) · B(Λc+ → K_S^0 X).")
  .with_decay_card(decay_card_signal)
  .apply
  .execute_on(scan_data + scan_incMC + exMC_signal)