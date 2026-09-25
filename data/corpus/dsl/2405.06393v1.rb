# Paper 2405.06393v1: e+e- -> p pbar pi0 at sqrt(s)=2.1000-3.0800 GeV
# Born cross section measurement using ConExc generator
# Partial reconstruction for low energies (2.1000-2.2324 GeV)
# Full reconstruction for higher energies (2.3094-3.0800 GeV)

# --- Datasets (BOSS 713 R-scan, 20 energy points) ---
data_2100  = DatasetManager.real_data.find("713_Rscan_2100")
data_2125  = DatasetManager.real_data.find("713_Rscan_2125")
data_2150  = DatasetManager.real_data.find("713_Rscan_2150")
data_2175  = DatasetManager.real_data.find("713_Rscan_2175")
data_2200  = DatasetManager.real_data.find("713_Rscan_2200")
data_2232  = DatasetManager.real_data.find("713_Rscan_2232")
data_2309  = DatasetManager.real_data.find("713_Rscan_2309")
data_2386  = DatasetManager.real_data.find("713_Rscan_2386")
data_2396  = DatasetManager.real_data.find("713_Rscan_2396")
data_2500  = DatasetManager.real_data.find("713_Rscan_2500")
data_2644  = DatasetManager.real_data.find("713_Rscan_2644")
data_2646  = DatasetManager.real_data.find("713_Rscan_2646")
data_2700  = DatasetManager.real_data.find("713_Rscan_2700")
data_2800  = DatasetManager.real_data.find("713_Rscan_2800")
data_2900  = DatasetManager.real_data.find("713_Rscan_2900")
data_2950  = DatasetManager.real_data.find("713_Rscan_2950")
data_2981  = DatasetManager.real_data.find("713_Rscan_2981")
data_3000  = DatasetManager.real_data.find("713_Rscan_3000")
data_3020  = DatasetManager.real_data.find("713_Rscan_3020")
data_3080  = DatasetManager.real_data.find("713_Rscan_3080")

data_points_low = [data_2100, data_2125, data_2150, data_2175, data_2200, data_2232]
data_points_high = [data_2309, data_2386, data_2396, data_2500, data_2644, data_2646,
                    data_2700, data_2800, data_2900, data_2950, data_2981,
                    data_3000, data_3020, data_3080]
data_points_all = data_points_low + data_points_high

# --- ConExc Decay Card (form D, -2 = user-supplied cross section) ---
# ConExc continuum e+e- -> p pbar pi0 with ISR and vacuum polarization
# PDG codes: 2212 = p+, -2212 = anti-p-, 111 = pi0
# xs_user.txt must be supplied externally with FDC-PWA amplitude model results
decay_card_ppbar_pi0 = <<~DECAYCARD
    Decay vpho
    1 ConExc -2 2212 -2212 111;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_ppbar_pi0 = DatasetManager.create_exclusive_mc_for(data_points_all) do |config|
  config.sample_name = "sig_ppbar_pi0_conexc"
  config.events = 200000
  config.decay_card = decay_card_ppbar_pi0
  config.cross_section = :default
end

# ==========================================================================
# Algorithm 1: Partial Reconstruction (Mode I) — missing pi0
# sqrt(s) = 2.1000 - 2.2324 GeV
# Reconstruct: 1 proton + 1 anti-proton, no photons
# Signal extracted from recoil mass M_recoil(p pbar) fit
# ==========================================================================
alg_mode_i = Algorithm.new("PPbarPi0PartialModeI")
alg_mode_i.set_header(["PPbarPi0PartialModeIAlg/PPbarPi0PartialModeI.h"])
           .set_constant({ "ECMS" => [:double, 2.1000] })
           .note(:partial_rec, "Mode I: missing pi0; reconstruct proton + anti-proton only")
           .note(:energy_range, "Applied at sqrt(s) = 2.1000-2.2324 GeV")
           .note(:theta_ppbar_cut, "Angle between proton and anti-proton < 173 degrees to suppress cosmic/beam background")
           .note(:cos_theta_recoil_cut, "|cos(theta_recoil)| <= 0.96 to suppress e+e- -> p pbar ISR background")
           .note(:fit_method, "Recoil mass M(p pbar)^recoil fitted with MC shape convolved with Gaussian + Landau background")
           .note(:conexc, "Born cross section measured using ConExc generator with ISR and VP corrections")

sel_i = Selection.new

sel_i.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
  nChrn ">=1"
}
.pid(method: :best) {
  identify :prp, against: [:pip, :kp]
  identify :prm, against: [:pim, :km]
}
.kinematic_fit([:prp, :prm]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_mode_i.with_decay_card(decay_card_ppbar_pi0).apply(sel_i)
alg_mode_i.execute_on(data_points_low + [exMC_ppbar_pi0])

# ==========================================================================
# Algorithm 2: Partial Reconstruction (Mode II) — missing proton
# sqrt(s) = 2.1000 - 2.2324 GeV
# Reconstruct: 1 anti-proton + at least 2 photons
# 1C kinematic fit for p_miss anti-p gamma gamma
# Signal from M(gamma gamma) pi0 mass fit
# ==========================================================================
alg_mode_ii = Algorithm.new("PPbarPi0PartialModeII")
alg_mode_ii.set_header(["PPbarPi0PartialModeIIAlg/PPbarPi0PartialModeII.h"])
             .set_constant({ "ECMS" => [:double, 2.1000] })
             .note(:missing_proton, "Mode II: missing proton; reconstruct anti-proton + 2 photons")
             .note(:energy_range, "Applied at sqrt(s) = 2.1000-2.2324 GeV")
             .note(:theta_gamma_ap_cut, "Angle between photon and closest anti-proton track > 30 degrees")
             .note(:kinematic_fit, "1C kinematic fit with missing proton, chi2_1C < 60")
             .note(:fit_method, "M(gamma gamma) fitted with MC shape convolved with Gaussian + linear background")

sel_ii = Selection.new

sel_ii.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrn ">=1"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 700
  angle_to_track 30.0
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=2"
}
.pid(method: :best) {
  identify :prm, against: [:pim, :km]
}
.kinematic_fit([:prm, :gamma, :gamma]) {
  nominal
  constrain_four_momentum
  chi2_cut 60
}

alg_mode_ii.with_decay_card(decay_card_ppbar_pi0).apply(sel_ii)
alg_mode_ii.execute_on(data_points_low + [exMC_ppbar_pi0])

# ==========================================================================
# Algorithm 3: Full Reconstruction
# sqrt(s) = 2.3094 - 3.0800 GeV
# Reconstruct: proton + anti-proton + at least 2 photons (pi0 -> gamma gamma)
# 4C kinematic fit, chi2_4C < 100
# Signal from M(gamma gamma) pi0 mass fit
# ==========================================================================
alg_full = Algorithm.new("PPbarPi0FullReco")
alg_full.set_header(["PPbarPi0FullRecoAlg/PPbarPi0FullReco.h"])
         .set_constant({ "ECMS" => [:double, 2.3094] })
         .note(:full_reconstruction, "Full reconstruction: proton + anti-proton + 2 photons")
         .note(:energy_range, "Applied at sqrt(s) = 2.3094-3.0800 GeV")
         .note(:theta_gamma_p_track, "Angle between photon and closest proton/anti-proton > 10/20 degrees")
         .note(:kinematic_fit, "4C kinematic fit for p anti-p gamma gamma; chi2_4C < 100; best combination by minimum chi2")
         .note(:fit_method, "M(gamma gamma) fitted with MC shape convolved with Gaussian + linear background")
         .note(:conexc, "Born cross section measured using ConExc generator with ISR and VP corrections")

sel_full = Selection.new

sel_full.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
  nChrn ">=1"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 700
  angle_to_track 10.0
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=2"
}
.pid(method: :best) {
  identify :prp, against: [:pip, :kp]
  identify :prm, against: [:pim, :km]
}
.kinematic_fit([:prp, :prm, :gamma, :gamma]) {
  nominal
  constrain_four_momentum
  chi2_cut 100
}

alg_full.with_decay_card(decay_card_ppbar_pi0).apply(sel_full)
alg_full.execute_on(data_points_high + [exMC_ppbar_pi0])