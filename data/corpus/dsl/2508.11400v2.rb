# 2508.11400v2: Lambda_c+ transverse polarization and decay asymmetry
# Single-tag approach, 5 decay modes, 13 energies 4600-4951 MeV, 6.4 fb^-1
# Uses DTagAlg single-tag (ST) pattern

DatasetManager.load_real_data("config/BES3_dataset.md")
DatasetManager.load_inclusive_mc("config/BES3_incMC.md")

scan_points = [
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4914"),
  DatasetManager.real_data.find("707_4946"),
]

incMC_points = scan_points.map { |d|
  DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
}

# Decay card: e+e- -> Lambda_c+ anti-Lambda_c- with inclusive anti-Lambda_c- decays
# Lambda_c+ decays via 5 channels; anti-Lambda_c- decays inclusively
decay_card = <<~DECAYCARD
  Decay psi(4520)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  0.200 p+ K_S0            PHSP;
  0.200 Lambda pi+         PHSP;
  0.200 Sigma0 pi+         PHSP;
  0.200 Sigma+ pi0         PHSP;
  0.200 p+ K- pi+          PHSP;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  Decay Lambda
  1.000 p+ pi- PHSP;
  Enddecay
  Decay Sigma0
  1.000 Lambda gamma PHSP;
  Enddecay
  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_Lambdac_pol"
  config.events        = 200_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

alg = TagAnalysis.new("LambdacPol")
alg.set_header(["LambdacPolAlg/LambdacPol.h"])
    .with_decay_card(decay_card)

# 5 tag modes: pK_S0, Lambda pi+, Sigma0 pi+, Sigma+ pi0, p K- pi+
# Sigma0 pi+ and Sigma+ pi0 have no vocabulary match -> dropped with note
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP,
          :LambdacPtoKsP,
          :LambdacPtoLambdaPi
end

alg.signal_side do |s|
  s.photons 0
end

alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg.note(:tag_mode_unavailable, "Lambda_c+ -> Sigma0 pi+ (Sigma0 -> Lambda gamma): no vocabulary match in DTagAlg enum. Omitted from tag modes.")
alg.note(:tag_mode_unavailable, "Lambda_c+ -> Sigma+ pi0 (Sigma+ -> p pi0): no vocabulary match in DTagAlg enum. Omitted from tag modes.")
alg.note(:multi_energy, "13 energy points 4600-4951 MeV. ECMS is placeholder; per-point CMS energy injected at runtime.")
alg.note(:polarization_analysis, "Transverse polarization and decay asymmetry parameters (alpha, beta, gamma) extracted from multidimensional angular analysis in ROOT on tagged Lambda_c+ samples.")
alg.note(:signal_yields, "Tag yields (from paper): pK-pi+ ~50083, pK_S0 ~9619, Lambda pi+ ~5742, Sigma0 pi+ ~2487, Sigma+ pi0 ~1268. Only first 3 modes available in DTagAlg.")

alg.apply
alg.execute_on(scan_points + incMC_points + sig_mc)