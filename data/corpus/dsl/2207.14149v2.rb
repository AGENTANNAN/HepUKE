# Paper: 2207.14149v2
# Title: Improved measurement of Lambda_c+ -> Lambda e+ nu_e and
#        first study of its internal dynamics
# Energy: 4.600-4.699 GeV (7 energy points)
# Double-tag (DT) method:
#   ST: Lambda_c- -> 14 hadronic tag modes
#   DT: Lambda_c+ -> Lambda e+ nu_e (semileptonic signal)
#   Lambda -> p pi-
# 4D ML fit for form factor extraction and |V_cs| determination

### Dataset preparation ###
data_703_4600 = DatasetManager.load_real_data.find("703_4600")

all_data = [data_703_4600]
all_incMC = DatasetManager.load_inclusive_mc

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Lc_to_Lambda_e_nu"
  config.events         = 500000
  config.decay_card     = nil  # TagAnalysis
end

### TagAnalysis: Double-tag with Lambda_c ST hadronic tags ###
alg = TagAnalysis.new("Lc_to_Lambda_e_nu")
alg.set_header(["LcLambdaENuAlg/LcLambdaENu.h"])

# ST: Lambda_c- -> 14 hadronic tag modes
alg.tag_side(:Lambdacm) {
  modes(:LambdacPtoPKPi, :LambdacPtoLambdaPi,
        :LambdacPtoLambdaPiPiPi, :LambdacPtoPKPiPiPi)
  charm -1
}

# Signal side: Lambda_c+ -> Lambda e+ nu_e (semileptonic)
# Lambda -> p pi-
# e+ identified via PID, missing neutrino
alg.signal_side {
  charged(ep: 1, prp: 1, pim: 1)
  missing :nu_e
}

# Fit: 4-momentum conservation with missing neutrino
alg.fit {
  constrain_four_momentum
  chi2_cut 200
}

# ST selection: M_BC and DeltaE cuts for Lambda_c- tags
# Total ST yield: 122,268 +/- 474
alg.note(:st_selection,
  "ST: 14 hadronic Lambda_c- tag modes. M_BC signal region fits. DeltaE requirements per mode. Total ST yield = 122,268+/-474. ST efficiencies from inclusive MC. Applied in ROOT.")

# DT selection: Lambda -> p pi-, e+ PID, U_miss variable for neutrino
# Background suppression: M(Lambda e+) < 2.12 GeV
alg.note(:dt_selection,
  "DT: Lambda_c+ -> Lambda e+ nu_e. Lambda: vertex-fit p pi-, |M-M(Lambda)| < 5 MeV. e+ PID: CL_e>0.001, CL_e>CL_pi, CL_e>CL_K. U_miss = E_miss - |p_miss|. U_miss fit: signal Gaussian+power-law tails, backgrounds from Lambda_c+->Lambda pi+pi0, Lambda_c+->Lambda mu+nu_mu, and combinatorial. DT yield = 1253+/-39. BF = (3.56+/-0.11+/-0.07)%. Applied in ROOT.")

# Form factor extraction: 4D ML fit to q2, cos_theta_e, cos_theta_p, chi
alg.note(:form_factors,
  "4D ML fit to kinematic variables (q2, cos_theta_e, cos_theta_p, chi) within -0.06 < U_miss < 0.06 GeV. Five free form-factor parameters: alpha1_gperp, alpha1_fperp, r_f+ = a0_f+/a0_gperp, r_fperp, r_g+. Determined a0_gperp = 0.54+/-0.04+/-0.01. First direct comparison of differential decay rates and form factors with LQCD. chi2/ndof = 0.85. Applied in ROOT.")

# |V_cs| extraction
alg.note(:vcs_extraction,
  "|V_cs| determination: combining measured BF with LQCD q2-integrated rate. |V_cs| = 0.936+/-0.017(BF)+/-0.024(LQCD)+/-0.007(tau_Lambdac). Consistent with |V_cs| from D->Klnu. Applied in ROOT.")

# 4.5 fb^-1 total at 4.600-4.699 GeV
alg.note(:energy_points,
  "7 c.m. energies: 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV. Total 4.5 fb^-1.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on([all_datasets])