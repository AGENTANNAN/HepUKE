# 2212.07214v2: Search for Lambda_c+ -> Sigma+ gamma at sqrt(s)=4.60-4.70 GeV
# TagAnalysis: Lambda_c- ST tags, signal side Lambda_c+ -> Sigma+ gamma
# Multi-energy: 7 energy points

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 7 energy points: 4.59953, 4.61184, 4.62800, 4.64067, 4.66122, 4.68184, 4.69857 GeV
# Using approximate sample names
energy_samples = %w[709_4600 709_4612 709_4628 709_4641 709_4661 709_4682 709_4699]

data_points = energy_samples.map { |s| DatasetManager.real_data.find(s) }
incMC_points = energy_samples.map { |s| DatasetManager.inclusive_mc.find(s) }

# Decay card: e+e- -> Lambda_c+ anti-Lambda_c-
# anti-Lambda_c- -> [tag mode] (ST)
# Lambda_c+ -> Sigma+ gamma, Sigma+ -> p pi0, pi0 -> gamma gamma
decay_card = <<~DECAY
Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- VSS;
Enddecay
Decay Lambda_c+
  1.000 anti-Sigma- gamma PHSP;
Enddecay
Decay anti-Sigma-
  1.000 anti-p- pi0 PHSP;
Enddecay
Decay anti-Lambda_c-
  1.000 anti-p- K+ pi- PHSP;
Enddecay
Decay pi0
  1.000 gamma gamma PHSP;
Enddecay
End
DECAY

exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_Lc_Sigma_gamma"
  config.events        = 500_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

alg = TagAnalysis.new("LcTagSigmaGamma", "00-00-01")
alg.set_header(["LcTagSigmaGamma/LcTagSigmaGamma.h"])
# Nominal ECMS — per-run MeasuredEcmsSvc handles energy variation for data
alg.set_constant({ "ECMS" => [:double, 4.600] })
alg.with_decay_card(decay_card)

# Tag side: anti-Lambda_c- reconstructed in hadronic modes (Table 2)
# Paper lists 10 tag modes:
#   pK-pi+, pKS0, Lambda pi-, pK+pi-pi0, pKS0pi0, Lambda pi-pi0,
#   pKS0 pi+pi-, Lambda pi- pi+pi-, Sigma0 pi-, Sigma- pi+pi-
# The three standard modes available in BOSS EvtRecDTag enum:
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP      # anti-Lambda_c- -> anti-p K+ pi-
  t.modes :LambdacPtoKsP        # anti-Lambda_c- -> anti-p KS0
  t.modes :LambdacPtoLambdaPi   # anti-Lambda_c- -> anti-Lambda pi-
end

alg.note(:tag_mode_unavailable,
  "Paper tag modes NOT in frozen BOSS EvtRecDTag enum: " \
  "pK+pi-pi0 (p K pi pi0), pKS0pi0, Lambda pi-pi0, pKS0 pi+pi-, " \
  "Lambda pi- pi+pi-, Sigma0 pi-, Sigma- pi+pi-. " \
  "These 7 of 10 modes are unavailable in BOSS 7.0.6/7.1.2. " \
  "Only 3 standard Lambda_c modes declared above. " \
  "Full mode list requires custom EvtRecDTag.h fork."
)

alg.note(:tag_side_selection,
  "ST selection: Charged tracks |cos(theta)|<0.93. Except for KS0/Lambda daughters, " \
  "|Vxy|<1cm, |Vz|<10cm. PID via dE/dx + TOF probability. " \
  "KS0 -> pi+pi-: |Vz|<20cm for daughters, secondary vertex fit (chi2<100), " \
  "decay length > 2x vertex resolution, M(pipi) in (0.487, 0.511) GeV/c^2. " \
  "Lambda -> anti-p pi+: PID on anti-p only, secondary vertex fit (chi2<100), " \
  "decay length > 2x vertex resolution, M(pbar pi+) in (1.111, 1.121) GeV/c^2. " \
  "Sigma0 -> gamma anti-Lambda: M in (1.179, 1.203) GeV/c^2. " \
  "Sigma- -> anti-p pi0: M in (1.176, 1.200) GeV/c^2. " \
  "pi0 -> gamma gamma: M in (0.115, 0.150) GeV/c^2, kinematic fit chi2<200. " \
  "DeltaE requirements per mode in Table 2. Candidates in M_BC signal region " \
  "(2.275, 2.310) GeV/c^2. Multiple candidates: keep minimum |DeltaE|."
)

# Signal side: Lambda_c+ -> Sigma+ gamma, Sigma+ -> p pi0
# One proton, one radiative gamma (>0.65 GeV), pi0 from gamma gamma
alg.signal_side do |s|
  s.photons 2                   # pi0 -> gamma gamma
  s.charged(prp: 1)             # exactly one proton from Sigma+ -> p pi0
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg.note(:signal_side_selection,
  "Signal side selection: exactly one additional charged track (proton), " \
  "|Vz|<20cm from IP. Proton from dE/dx+TOF probability. " \
  "Radiative photon energy > 0.65 GeV. Exactly one good photon for signal. " \
  "pi0 photon energies < 0.45 GeV. pi0 from gamma gamma: " \
  "M(gammagamma) in (0.115, 0.150) GeV/c^2. " \
  "Sigma+ candidate: M(p pi0) in (1.176, 1.200) GeV/c^2. " \
  "DeltaE_sig in (-0.038, 0.026) GeV. " \
  "Multiple pi0 candidates: keep minimum |DeltaE_sig|."
)

# Fit: 4-momentum conservation
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply
alg.execute_on(data_points + incMC_points + exMCs)