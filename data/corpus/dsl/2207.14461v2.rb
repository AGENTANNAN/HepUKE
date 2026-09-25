# Paper: 2207.14461v2
# Title: Measurement of the absolute branching fraction of the singly
#        Cabibbo suppressed decay Lambda_c+ -> p eta'
# Energy: 4.600-4.699 GeV (7 energy points)
# Double-tag (DT) method:
#   ST: Lambda_c- -> 10 hadronic tag modes
#   DT: Lambda_c+ -> p eta'
#   eta' -> pi+ pi- eta (eta -> gamma gamma) and eta' -> pi+ pi- gamma
# 4.5 fb^-1 total

### Dataset preparation ###
data_703_4600 = DatasetManager.load_real_data.find("703_4600")

all_data = [data_703_4600]
all_incMC = DatasetManager.load_inclusive_mc

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Lc_to_p_etap"
  config.events         = 500000
  config.decay_card     = nil  # TagAnalysis
end

### TagAnalysis: Double-tag with Lambda_c ST hadronic tags ###
alg = TagAnalysis.new("Lc_to_p_etap")
alg.set_header(["LcpEtapAlg/LcpEtap.h"])

# ST: Lambda_c- -> 10 hadronic tag modes
# Note: paper uses 10 tag modes (including Sigma modes). The 4 family-level
# symbols below cover the main decay topologies; finer-grained modes handled in ROOT.
alg.tag_side(:Lambdacm) {
  modes(:LambdacPtoPKPi, :LambdacPtoLambdaPi,
        :LambdacPtoLambdaPiPiPi, :LambdacPtoPKPiPiPi)
  charm -1
}

# Signal side: Lambda_c+ -> p eta'
# eta' -> pi+ pi- eta (eta->gamma gamma) or eta' -> pi+ pi- gamma
alg.signal_side {
  charged(prp: 1, pip: 1, pim: 1)
  photons 2
}

# Fit
alg.fit {
  constrain_four_momentum
  chi2_cut 200
}

# ST selection
alg.note(:st_selection,
  "ST: 10 hadronic Lambda_c- tag modes. M_BC signal region (2.275,2.310) GeV/c2. Asymmetric DeltaE requirements per mode. ST yields from unbinned ML fits to M_BC (signal: MC shape convolved with Gaussian; background: ARGUS). Sum of ST yields = (1.0524+/-0.0038)e5. Applied in ROOT.")

# DT selection for eta' -> pi+ pi- eta and eta' -> pi+ pi- gamma
alg.note(:dt_selection,
  "DT: Lambda_c+ -> p eta'. Exactly 3 tight charged tracks recoiling against ST Lambda_c-. PID: proton, pi+, pi-. Lambda veto: M(p pi-)>1.125 GeV/c2. K_S0 veto: M(pi+pi-) not in (0.490,0.510) GeV/c2. eta' -> pi+ pi- eta: eta treated as missing, recoil mass Mrec(Lambdac- p pi+ pi-) in eta window. eta' -> pi+ pi- gamma: gamma from unused photons, DeltaE selection. Simultaneous unbinned ML fit to Mrec(Lambdac- p) and M(pi+pi-gamma) for BF extraction. N_sig = 9.2 events total. BF = (5.62+2.46/-2.04+/-0.26)e-4. Significance 3.6sigma. Applied in ROOT.")

# 7 energy points
alg.note(:energy_points,
  "7 c.m. energies: 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV. Total 4.5 fb^-1.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on([all_datasets])