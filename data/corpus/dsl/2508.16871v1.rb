# 2508.16871v1: e+e- -> Sigma_c anti-Sigma_c and Lambda_c+ anti-Sigma_c- upper limits
# Tag method with 3 sample types: Tag Lambda_c+, Tag Lambda_c+ pi, Tag Lambda_c+ pi pi
# 4 Lambda_c+ decay modes: pK_S0, pK-pi+, Lambda pi+, Sigma0 pi+
# 2 energies for Sigma_c Sigma_c (4.918, 4.951 GeV)
# 5 energies for Lambda_c+ anti-Sigma_c- (4.750, 4.781, 4.843, 4.918, 4.951 GeV)

DatasetManager.load_real_data("config/BES3_dataset.md")
DatasetManager.load_inclusive_mc("config/BES3_incMC.md")

scan_points_sigmac_sigmac = [
  DatasetManager.real_data.find("707_4914"),
  DatasetManager.real_data.find("707_4946"),
]

scan_points_lc_sigmac = [
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4914"),
  DatasetManager.real_data.find("707_4946"),
]

all_points = (scan_points_sigmac_sigmac + scan_points_lc_sigmac).uniq { |d| "#{d.boss}_#{d.sample_name}" }

incMC_points = all_points.map { |d|
  DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
}

decay_card = <<~DECAYCARD
  Decay psi(4520)
  1.000 Lambda_c+ anti-Sigma_c- PHSP;
  Enddecay
  Decay Lambda_c+
  0.250 p+ K_S0            PHSP;
  0.250 p+ K- pi+          PHSP;
  0.250 Lambda pi+         PHSP;
  0.250 Sigma0 pi+         PHSP;
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
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc_for(all_points) do |config|
  config.sample_name   = "sig_Lc_Sigmac"
  config.events        = 200_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

# ---- a) Tag Lambda_c+ sample: single tag ST, search Sigma_c in signal side ----

alg_lc_st = TagAnalysis.new("LambdacSigmaC_ST")
alg_lc_st.set_header(["LambdacSigmaC_STAlg/LambdacSigmaC_ST.h"])
    .with_decay_card(decay_card)

alg_lc_st.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP,
          :LambdacPtoKsP,
          :LambdacPtoLambdaPi
end

# Signal side: no tagged charged tracks at BOSS level.
# Sigma_c reconstructed via additional pion(s) from stored tag-side variables in ROOT
# RM(Lambda_c+)+M(Lambda_c+)-m(Lambda_c+) > 2.54 GeV, M_BC(Lambda_c+pi) in ROOT
alg_lc_st.signal_side do |s|
  s.photons 0
end

alg_lc_st.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg_lc_st.note(:tag_mode_unavailable, "Lambda_c+ -> Sigma0 pi+ (Sigma0 -> Lambda gamma): no vocabulary match in DTagAlg enum. Omitted from tag modes.")
alg_lc_st.note(:sigma_c_reconstruction, "Sigma_c(2455/2520) reconstructed in ROOT via M_BC(Lambda_c+ pi) and RM(Lambda_c+) > 2.54 GeV/c^2 after requiring |DeltaM| < 0.02 GeV/c^2 and Lambda_c+ vertex fit chi2 < 200. Not expressible in BOSS DSL.")
alg_lc_st.note(:three_sample_types, "Three tag samples: Tag Lambda_c+ (ST), Tag Lambda_c+ pi (ST+1pi), Tag Lambda_c+ pi pi (ST+2pi). Only the ST sample is expressed in DSL; the additional pion combinations are analyzed in ROOT.")
alg_lc_st.note(:upper_limits, "No significant signal observed. Upper limits at 90% C.L. set in ROOT analysis via simultaneous fit.")
alg_lc_st.note(:dataset_validation, "7 energy points total (2 for Sigma_c Sigma_c, 5 for Lambda_c+ anti-Sigma_c-). Exact energy values validated against paper.")

alg_lc_st.apply
alg_lc_st.execute_on(all_points + incMC_points + sig_mc)