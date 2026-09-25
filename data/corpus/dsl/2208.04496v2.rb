# Paper: 2208.04496v2
# Title: Search for a massless dark photon gamma' in Lambda_c+ -> p gamma'
# Energy: 4.600-4.699 GeV (7 energy points)
# Double-tag (DT) method:
#   ST: Lambda_c- -> 10 hadronic tag modes
#   DT: Lambda_c+ -> p gamma' (gamma' invisible, missing energy)
# First search for FCNC dark photon in charmed baryon sector
# 4.5 fb^-1 total

### Dataset preparation ###
data_703_4600 = DatasetManager.load_real_data.find("703_4600")

all_data = [data_703_4600]
all_incMC = DatasetManager.load_inclusive_mc

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Lc_to_p_gammap"
  config.events         = 500000
  config.decay_card     = nil  # TagAnalysis
end

### TagAnalysis: Double-tag with Lambda_c ST hadronic tags ###
alg = TagAnalysis.new("Lc_to_p_gammap")
alg.set_header(["LcpGammapAlg/LcpGammap.h"])

# ST: Lambda_c- -> 10 hadronic tag modes
# Note: paper uses 10 tag modes. The 4 family-level symbols below
# cover the main decay topologies; finer-grained modes handled in ROOT.
alg.tag_side(:Lambdacm) {
  modes(:LambdacPtoPKPi, :LambdacPtoLambdaPi,
        :LambdacPtoLambdaPiPiPi, :LambdacPtoPKPiPiPi)
  charm -1
}

# Signal side: Lambda_c+ -> p gamma'
# Proton from Lambda_c+ decay; gamma' is massless and invisible
# Modeled as missing :nu (generic massless invisible particle)
alg.signal_side {
  charged(prp: 1)
  missing :nu
}

# Fit
alg.fit {
  constrain_four_momentum
  chi2_cut 200
}

# ST selection
alg.note(:st_selection,
  "ST: 10 hadronic Lambda_c- tag modes. M_BC signal region (2.275,2.310) GeV/c2. Asymmetric DeltaE requirements per mode. ST yields from unbinned ML fits to M_BC (signal: MC shape convolved with Gaussian; background: ARGUS). Sum of ST yields = 105,244+/-384. Applied in ROOT.")

# DT selection: Lambda_c+ -> p gamma'
alg.note(:dt_selection,
  "DT: Lambda_c+ -> p gamma'. Only 1 tight track recoiling against ST Lambda_c-, identified as proton. No loose tracks allowed. gamma' is invisible: Emax < 0.3 GeV and Esum < 0.5 GeV on unused showers to veto pi0 background. Signal region: M^2_rec(Lambdac- p) in (0.0, 0.1) GeV^2/c^4. Main background: Lambda_c+ -> p K_L0 (peaking) and Lambda_c+ Lambda_c- backgrounds. N_obs = 13, N_bkg = 14.6+/-1.5. Upper limit B(Lambda_c+ -> p gamma') < 8.0e-5 at 90% CL. Applied in ROOT.")

# 7 energy points
alg.note(:energy_points,
  "7 c.m. energies: 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV. Total 4.5 fb^-1.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on([all_datasets])