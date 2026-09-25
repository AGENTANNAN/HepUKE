# Paper: 2206.10791v1
# Title: Measurement of Lambda -> n gamma decay branching fraction and asymmetry
# Energy: 3.097 GeV (J/psi, single energy point)
# Final state: J/psi -> Lambda Lambdabar; Lambdabar -> pbar pi+ (ST); Lambda -> n gamma (DT, neutron missing)
# Uses double-tag technique: reconstruct ST Lambdabar first, then signal Lambda -> n gamma
# Two charge-conjugate modes: Lambda -> n gamma and Lambdabar -> nbar gamma

### Dataset preparation ###
data_708_3097 = DatasetManager.load_real_data.find("708_3097")
incMC_708_3097 = DatasetManager.load_inclusive_mc.find("708_3097")

all_data = [data_708_3097]
all_incMC = [incMC_708_3097]

# Decay card for ST: J/psi -> Lambda Lambdabar, Lambdabar -> pbar pi+
# and signal side: Lambda -> n gamma
decay_card_st = <<~DECAYCARD
    Decay J/psi
    1.0000  anti-Lambda0  Lambda0     HELAMP 1.0 0.0 0.461 0.0;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+                HypWK;
    Enddecay

    Decay Lambda0
    1.0000  n0  gamma                   PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Jpsi_LL_ngamma_signal"
  config.events         = 200000
  config.decay_card     = decay_card_st
  config.cross_section  = :default
end

### Event selection ###
# ST side: reconstruct Lambdabar -> pbar pi+
# Signal side: Lambda -> n gamma with neutron as missing particle
# Two modes: mode A (Lambda -> n gamma, 1C fit with n missing)
#            mode B (Lambdabar -> nbar gamma, 3C equivalent, nbar direction measured)

alg_mode_A = Algorithm.new("Lambda_nGamma_1C")
alg_mode_A.set_header(["Lambda_nGamma_1CAlg/Lambda_nGamma_1C.h"])
           .set_alias({"std::vector<double>" => "Vdouble"})

sel_mode_A = Selection.new

sel_mode_A.select_track {
                 cos_theta 0.93
                 nTot ">=2"    # at least pbar + pi+ from ST
               }
               .select_photon {
                 tdc_emc_start 0
                 tdc_emc_end   700
                 angle_to_track 10.0
                 energyThreshold_b 0.025
                 energyThreshold_e 0.050
                 nGam ">=1"
               }
               # PID: anti-proton identification via largest likelihood
               # Pions from Lambdabar have p < 0.5 GeV/c, anti-protons p > 0.5 GeV/c
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:kaon, :pion]
                 nprm ">=1"    # anti-proton
               }
               .remove([:prm <= :chrgn])
               # Remaining tracks assigned as pions
               .assign({:chrgn => :pip})
               # Select isolated photon for signal side (photon from Lambda -> n gamma)
               # Anti-proton angular isolation for photons (20 deg for anti-proton tracks)
               .select_isolated_photon {
                 angle_to_prm_track 20.0
                 nGam ">=1"
               }
               # ST: reconstruct Lambdabar via secondary vertex fit of pbar pi+
               .secondary_vertex_fit([:prm, :pip]) {
                 build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                 chi2_cut 20
                 remove_used_particle_from_candidate_list
               }

alg_mode_A.with_decay_card(decay_card_st).apply(sel_mode_A)

# 1C kinematic fit: J/psi -> Lambdabar(Lambda_bar) Lambda n gamma
# With neutron as missing particle
# Note: this is applied in the analysis at the signal side level,
# after ST Lambdabar reconstruction.
alg_mode_A.note(:kinematic_fit_1C,
  "1C kinematic fit: J/psi -> Lambda_bar(Lambda_bar) Lambda n gamma with neutron missing. chi2_1C < 10. Applied in ROOT.")

# Pi0 veto: 1C/3C kinematic fit under J/psi -> Lambda_bar n gamma gamma hypothesis,
# reject events with M(gammagamma) within 20 MeV/c^2 of pi0 mass
alg_mode_A.note(:pi0_veto,
  "Pi0 veto via kinematic fit under J/psi -> Lambda_bar n gamma gamma. M(gammagamma) within 20 MeV/c^2 of pi0 mass rejected. Applied in ROOT.")

# Lambda/Lambdabar mass window: 8 MeV/c^2 around nominal mass
alg_mode_A.note(:lambda_mass_window,
  "M(pbar pi+) within 8 MeV/c^2 of Lambda mass; M_rec(Lambda_bar) within [1.03, 1.18] GeV/c^2. Applied in ROOT.")

# Decay length > 2*sigma for Lambda candidates
alg_mode_A.note(:decay_length,
  "Decay length > 2*vertex_resolution for Lambda candidates; applied in ROOT.")

# BDT on photon to discriminate signal from background
alg_mode_A.note(:bdt_photon,
  "BDT on photon (deposited energy, secondary moment, nhits, Zernike A42, shape). BDT > 0.3 required. Applied in ROOT.")

# Photon energy > 150 MeV, opening angle > 20 deg from neutron candidate (BG B suppression)
alg_mode_A.note(:photon_quality_cuts,
  "Photon E > 150 MeV; opening angle between photon and neutron > 20 deg. Applied in ROOT.")

# Photon E < 400 MeV for signal photon
alg_mode_A.note(:photon_energy_upper,
  "Signal photon E < 400 MeV; most energetic shower > 0.4 GeV = neutron candidate. Applied in ROOT.")

# Signal extracted via fit to E_gamma_Lambda distribution in ROOT
alg_mode_A.note(:signal_extraction,
  "Signal yield from unbinned max likelihood fit to E_gamma^{Lambda} distribution. BG modeled with MC shapes convolved with Gaussian. Applied in ROOT.")

all_datasets = all_data + all_incMC + exMC_signal
root_files_A = alg_mode_A.execute_on(all_datasets)