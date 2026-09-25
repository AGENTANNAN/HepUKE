### Dataset preparation ###
# Seven BESIII energy points between 4599.53 MeV and 4698.82 MeV (total ~4.5 fb^-1)
data_4600 = DatasetManager.real_data.find("703_4600")   # 4599.53 MeV  (~587 pb^-1)
data_4610 = DatasetManager.real_data.find("706_4610")   # 4611.86 MeV  (~104 pb^-1)
data_4620 = DatasetManager.real_data.find("706_4620")   # 4628.00 MeV  (~522 pb^-1)
data_4640 = DatasetManager.real_data.find("706_4640")   # 4640.91 MeV  (~552 pb^-1)
data_4660 = DatasetManager.real_data.find("706_4660")   # 4661.24 MeV  (~529 pb^-1)
data_4680 = DatasetManager.real_data.find("706_4680")   # 4681.92 MeV  (~1667 pb^-1)
data_4700 = DatasetManager.real_data.find("706_4700")   # 4698.82 MeV  (~536 pb^-1)

# Corresponding inclusive MC samples
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

data_points  = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]
incMC_points = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# --- Phase-space exclusive signal MC : Lambda_c+ -> n K_S^0 pi+ ---
decay_card_nKsPi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 n0 K_S0 pi+ PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-n0 K_S0 pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# --- Phase-space exclusive signal MC : Lambda_c+ -> n K_S^0 K+ ---
decay_card_nKsK = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 n0 K_S0 K+ PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-n0 K_S0 K- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Signal MC generated at every energy point (same mode, different related dataset)
exMC_nKsPi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_nKsPi"
  config.events        = 100000
  config.decay_card    = decay_card_nKsPi
  config.cross_section = :default
end

exMC_nKsK = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_nKsK"
  config.events        = 100000
  config.decay_card    = decay_card_nKsK
  config.cross_section = :default
end

### Tag-based event selection (Lambda_c double tag) ###
# Ten hadronic tag modes (the eleventh mode is omitted: not in the frozen tag-mode list)
tag_modes = [
  :LambdacPtoPKs,          # p K_S^0
  :LambdacPtoPKPi,         # p K+ pi-
  :LambdacPtoPKsPi0,       # p K_S^0 pi0
  :LambdacPtoPKsPiPi,      # p K_S^0 pi+ pi-
  :LambdacPtoPKPiPi0,      # p K+ pi- pi0
  :LambdacPtoLambdaPi,     # Lambda pi-
  :LambdacPtoLambdaPiPi0,  # Lambda pi- pi0
  :LambdacPtoLambdaPiPiPi, # Lambda pi- pi+ pi-
  :LambdacPtoSigma0Pi,     # Sigma0 pi+  (Sigma0 -> Lambda gamma)
  :LambdacPtoSigmamPiPi    # Sigma- pi+ pi+
]

# ---------------- Mode I : Lambda_c+ -> n K_S^0 pi+ ----------------
alg_nKsPi = TagAnalysis.new("LcDTnKsPi")
alg_nKsPi.set_header(["LcDTnKsPiAlg/LcDTnKsPi.h"])
         .set_constant({ "ECMS" => [:double, 4.640] })  # representative energy of the 7-point set
         .with_decay_card(decay_card_nKsPi)

# Tag side: Lambda_c (charm = -1) in the ten hadronic modes
alg_nKsPi.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm(-1)
end

# Signal side: K_S^0 -> pi+pi- (2 charged pions) + bachelor pi+, and a missing neutron
alg_nKsPi.signal_side do |s|
  s.charged(pip: 2, pim: 1)
  s.missing :n                        # neutron identified via Mmiss^2 = E_miss^2 - p_miss^2
  s.min_photon_angle 10.0
end

# Kinematic fit: 4-momentum constraint, chi2 < 200
alg_nKsPi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_nKsPi
  .note(:ks0_secondary_vertex, "K_S^0 reconstructed from the two signal-side pions via a secondary vertex fit; require flight-length significance L/sigma_L > 2 and keep the candidate with the largest L/sigma_L")
  .note(:loose_daughter_pid, "K_S^0 daughters kept under a loose PID requirement chi2_pi < 4 OR chi2_K < 4")
  .note(:bachelor_selection, "bachelor pi+ required with |Vz| < 10 cm and Vxy < 1 cm and a loose PID requirement; likelihood separation L(pi) > L(K) and L(pi) > 0")
  .note(:missing_neutron, "missing neutron identified through Mmiss^2 = E_miss^2 - p_miss^2, expected near the neutron mass squared")
  .note(:background_veto, "veto 0.235 < M_npi+ - Mmiss < 0.265 GeV/c^2 and 0.240 < M_npi- - Mmiss < 0.270 GeV/c^2 (mode n K_S^0 pi+)")
  .note(:omitted_mode, "the eleventh Lambda_c tag mode is omitted because it is not present in the frozen tag-mode list")
  .note(:yield_extraction, "yields extracted from a 2D unbinned maximum-likelihood fit to Mmiss^2 versus M(pi+pi-) using PHSP MC signal shapes and ST-sideband plus peaking Lambda_c+ background descriptions")

alg_nKsPi.apply
root_files_nKsPi = alg_nKsPi.execute_on(data_points + incMC_points + exMC_nKsPi)

# ---------------- Mode II : Lambda_c+ -> n K_S^0 K+ ----------------
alg_nKsK = TagAnalysis.new("LcDTnKsK")
alg_nKsK.set_header(["LcDTnKsKAlg/LcDTnKsK.h"])
        .set_constant({ "ECMS" => [:double, 4.640] })
        .with_decay_card(decay_card_nKsK)

# Tag side identical to Mode I
alg_nKsK.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm(-1)
end

# Signal side: K_S^0 -> pi+pi- (2 charged pions) + bachelor K+, and a missing neutron
alg_nKsK.signal_side do |s|
  s.charged(pip: 1, pim: 1, kp: 1)
  s.missing :n                        # neutron identified via Mmiss^2 = E_miss^2 - p_miss^2
  s.min_photon_angle 10.0
end

# Kinematic fit: 4-momentum constraint, chi2 < 200
alg_nKsK.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_nKsK
  .note(:ks0_secondary_vertex, "K_S^0 reconstructed from the two signal-side pions via a secondary vertex fit; require flight-length significance L/sigma_L > 2 and keep the candidate with the largest L/sigma_L")
  .note(:loose_daughter_pid, "K_S^0 daughters kept under a loose PID requirement chi2_pi < 4 OR chi2_K < 4")
  .note(:bachelor_selection, "bachelor K+ required with |Vz| < 10 cm and Vxy < 1 cm and a loose PID requirement; likelihood separation L(K) > L(pi) and L(K) > 0")
  .note(:missing_neutron, "missing neutron identified through Mmiss^2 = E_miss^2 - p_miss^2, expected near the neutron mass squared")
  .note(:background_veto, "veto 0.240 < M_npi+ - Mmiss < 0.260 GeV/c^2 and 0.248 < M_npi- - Mmiss < 0.268 GeV/c^2 (mode n K_S^0 K+)")
  .note(:omitted_mode, "the eleventh Lambda_c tag mode is omitted because it is not present in the frozen tag-mode list")
  .note(:yield_extraction, "yields extracted from a 2D unbinned maximum-likelihood fit to Mmiss^2 versus M(pi+pi-) using PHSP MC signal shapes and ST-sideband plus peaking Lambda_c+ background descriptions")

alg_nKsK.apply
root_files_nKsK = alg_nKsK.execute_on(data_points + incMC_points + exMC_nKsK)