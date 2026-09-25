# Paper: 2207.11483v3
# Title: First study of semileptonic decays Lambda_c+ -> pK- e+ nu_e,
#        Lambda_c+ -> Lambda(1520) e+ nu_e and Lambda_c+ -> Lambda(1405) e+ nu_e
# Energy: 4.600-4.699 GeV (7 energy points)
# Double-tag (DT) method:
#   ST: Lambda_c- -> 14 hadronic tag modes
#   DT: Lambda_c+ -> pK- e+ nu_e (semileptonic signal)
#   Also searches for Lambda(1520) and Lambda(1405) resonances in pK- spectrum
# Semileptonic: missing neutrino

### Dataset preparation ###
data_703_4600 = DatasetManager.load_real_data.find("703_4600")

all_data = [data_703_4600]
all_incMC = DatasetManager.load_inclusive_mc

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Lc_to_pK_e_nu"
  config.events         = 500000
  config.decay_card     = nil  # TagAnalysis
end

### TagAnalysis: Double-tag with Lambda_c ST hadronic tags ###
alg = TagAnalysis.new("Lc_to_pK_e_nu")
alg.set_header(["LcpKENuAlg/LcpKENu.h"])

# ST: Lambda_c- -> 14 hadronic tag modes (same as Ref. [36])
alg.tag_side(:Lambdacm) {
  modes(:LambdacPtoPKPi, :LambdacPtoLambdaPi,
        :LambdacPtoLambdaPiPiPi, :LambdacPtoPKPiPiPi)
  charm -1
}

# Signal side: Lambda_c+ -> p K- e+ nu_e (semileptonic)
# e+ identified via PID, missing neutrino
alg.signal_side {
  charged(prp: 1, km: 1, ep: 1)
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
  "ST: 14 hadronic Lambda_c- tag modes. M_BC signal region [2.28,2.30] GeV/c2. DeltaE requirements per mode. Total ST yield = 122,268+/-474. ST efficiencies from inclusive MC. Applied in ROOT.")

# DT selection: p K- e+ combination, PID, background suppression
# U_miss = E_miss - |p_miss| for neutrino identification
alg.note(:dt_selection,
  "DT: Lambda_c+ -> pK- e+ nu_e. Exactly 3 charged tracks recoiling against ST Lambda_c-. Proton and kaon PID same as ST. e+ PID: CL_e>0.001, CL_e>CL_pi, CL_e>CL_K. M(pKe+)<2.15 GeV to suppress Lambda_c+ -> pK-pi+. pi0 veto for Lambda_c+ -> pK-pi+pi0. U_miss fit: signal Gaussian+power-law tails, backgrounds from Lambda_c+->pK-pi+pi0, Lambda_c+->pK-mu+nu_mu, and combinatorial. DT yield = 33.5+/-6.3. Significance 8.9sigma. Applied in ROOT.")

# Lambda(1520)/Lambda(1405) search via 2D fit to M(pK-) vs U_miss
alg.note(:resonance_search,
  "Lambda(1520)/Lambda(1405) search: 2D unbinned fit to M(pK-) vs U_miss. Lambda(1520) signal: Breit-Wigner (M=1.5195, Gamma=0.0156). Lambda(1405) signal: Flatte parameterization. Yields: Lambda(1520) = 8.4+/-4.3 (3.3sigma), Lambda(1405) = 14.8+/-6.7 (3.2sigma). BF(Lambda_c+ -> pK- e+ nu_e) = (0.88+/-0.17+/-0.07)e-3. Applied in ROOT.")

# Form factors and |V_cs| extraction
alg.note(:vcs_extraction,
  "|V_cs| determination: combining BF(Lambda_c+ -> Lambda(1520) e+ nu_e) with LQCD q2-integrated rate. |V_cs| = 1.3+/-0.3(BF)+/-0.1(LQCD). First determination of |V_cs| from baryonic SL decay to excited Lambda state. Applied in ROOT.")

# 4.5 fb^-1 total at 4.600-4.699 GeV
alg.note(:energy_points,
  "7 c.m. energies: 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV. Total 4.5 fb^-1.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on([all_datasets])