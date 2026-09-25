# Paper 2404.13840v2: Study of e+e- -> omega X(3872) and e+e- -> gamma X(3872) at 4.66-4.95 GeV
# Three analysis sections: A (omega gamma J/psi), B (omega pi+pi-J/psi), C (gamma pi+pi-J/psi)
# Uses lepton momentum/EMC based e/mu separation

# --- Datasets (BOSS 706 + 707, 4.66-4.95 GeV) ---
# Energy points covering 4.66-4.95 GeV range, total ~4.5 fb^-1
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")

data_points = [data_4660, data_4680, data_4700,
               data_4740, data_4750, data_4780, data_4840]

incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_4740 = DatasetManager.inclusive_mc.find("707_4740")
incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")

incMC_points = [incMC_4660, incMC_4680, incMC_4700,
                incMC_4740, incMC_4750, incMC_4780, incMC_4840]

# ==========================================================================
# Section A: omega gamma J/psi events
# Final state: pi+ pi- l+ l- + at least 2 out of 3 photons (1 photon undetected)
# omega -> pi+ pi- pi0 (pi0 -> gamma gamma)
# X(3872) -> gamma J/psi OR chi_cJ -> gamma J/psi
# J/psi -> l+ l- (l = e, mu)
# ==========================================================================

# Decay card for omega X(3872) -> omega gamma J/psi channel
decay_card_omega_gamma_jpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000 omega X(3872) PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay X(3872)
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for omega chi_cJ -> omega gamma J/psi channel
decay_card_omega_chi_cj_gamma_jpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000 omega chi_c1 PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

exMC_omega_gamma_jpsi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_omegaX3872_gammaJpsi"
  config.events = 200000
  config.decay_card = decay_card_omega_gamma_jpsi
  config.cross_section = :default
end

exMC_omega_chi_cj = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_omegaChiCJ_gammaJpsi"
  config.events = 200000
  config.decay_card = decay_card_omega_chi_cj_gamma_jpsi
  config.cross_section = :default
end

# --- Algorithm A: omega gamma J/psi ---
alg_a = Algorithm.new("OmegaGammaJpsi")
alg_a.set_header(["OmegaGammaJpsiAlg/OmegaGammaJpsi.h"])
     .set_constant({ "ECMS" => [:double, 4.660] })
     .note(:multi_energy, "Data taken at multiple energy points: 4.61-4.95 GeV; run separately per dataset")
     .note(:lepton_separation, "Charged tracks with p<1 GeV/c assigned as pions, p>1 GeV/c as leptons; EMC deposit <0.4 GeV for muons, >1.0 GeV for electrons")
     .note(:missing_photon, "One of three photons is undetected; 1C kinematic fit with missing mass constrained to zero")
     .note(:two_c_fit, "2C kinematic fit constraining pi0 mass for photon assignment; chi2<20 selected")
     .note(:psi3686_veto, "Veto RM(pi+pi-) within 0.03 GeV of psi(3686) mass and |M(pi+pi-J/psi)-m_psi3686|>0.01 GeV")
     .note(:muon_muc_depth, "At least one muon candidate must have MUC hit depth >30 cm to suppress mu/pi mis-ID")
     .note(:omega_signal, "M(pi+pi-pi0) in [0.74,0.82] GeV/c^2; sideband [0.64,0.72] U [0.84,0.92]")
     .note(:jpsi_signal, "M(l+l-) in [3.06,3.14] GeV/c^2; sideband [2.96,3.04] U [3.17,3.25]")

sel_a = Selection.new

sel_a.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 700
  angle_to_track 10.0
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=2"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :pion, against: [:kaon]
  identify :electron, against: [:pion]
  identify :muon, against: [:pion]
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=0"
}
.kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :ep, :em]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_a.with_decay_card(decay_card_omega_gamma_jpsi).apply(sel_a)
alg_a.execute_on(data_points + incMC_points + [exMC_omega_gamma_jpsi, exMC_omega_chi_cj])

# ==========================================================================
# Section B: omega pi+ pi- J/psi events (X(3872)->pi+pi-J/psi)
# 5-track mode: pi+ pi- pi+- l+ l- gamma gamma
# 6-track mode: pi+ pi- pi+ pi- l+ l- gamma (at least 1 photon)
# ==========================================================================

decay_card_omega_pipi_jpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000 omega X(3872) PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay X(3872)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

exMC_omega_pipi_jpsi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_omegaX3872_pipiJpsi"
  config.events = 200000
  config.decay_card = decay_card_omega_pipi_jpsi
  config.cross_section = :default
end

alg_b = Algorithm.new("OmegaPiPiJpsi")
alg_b.set_header(["OmegaPiPiJpsiAlg/OmegaPiPiJpsi.h"])
     .set_constant({ "ECMS" => [:double, 4.660] })
     .note(:multi_energy, "Data taken at multiple energy points: 4.61-4.95 GeV")
     .note(:lepton_separation, "Charged tracks with p<1 GeV/c assigned as pions, p>1 GeV/c as leptons; EMC <0.4 GeV for muons, >1.0 GeV for electrons")
     .note(:five_track_mode, "5-track events: pi+ pi- pi+- l+ l- + 2 photons; pi0 from gamma gamma")
     .note(:six_track_mode, "6-track events: pi+ pi- pi+ pi- l+ l- + at least 1 photon; pi+ pi- from X(3872), pi+ pi- + pi0 from omega")
     .note(:omega_signal, "M(pi+pi-pi0) in [0.75,0.81] GeV/c^2; sideband [0.66,0.72] U [0.84,0.90]")
     .note(:jpsi_signal, "M(l+l-) in [3.07,3.13] GeV/c^2; sideband [2.98,3.04] U [3.17,3.23]")
     .note(:r_measurement, "R = B(X(3872)->gamma J/psi) / B(X(3872)->pi+pi-J/psi) measured from omega X(3872) events")

sel_b = Selection.new

sel_b.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 700
  angle_to_track 10.0
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=1"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :pion, against: [:kaon]
  identify :electron, against: [:pion]
  identify :muon, against: [:pion]
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=0"
}
.kinematic_fit([:pip, :pim, :pip, :pim, :ep, :em, :gamma, :gamma]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_b.with_decay_card(decay_card_omega_pipi_jpsi).apply(sel_b)
alg_b.execute_on(data_points + incMC_points + [exMC_omega_pipi_jpsi])

# ==========================================================================
# Section C: gamma pi+ pi- J/psi events (gamma X(3872)->pi+pi-J/psi)
# Final state: pi+ pi- l+ l- + at least 1 photon
# ==========================================================================

decay_card_gamma_pipi_jpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000 gamma X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

exMC_gamma_pipi_jpsi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_gammaX3872_pipiJpsi"
  config.events = 200000
  config.decay_card = decay_card_gamma_pipi_jpsi
  config.cross_section = :default
end

alg_c = Algorithm.new("GammaPiPiJpsi")
alg_c.set_header(["GammaPiPiJpsiAlg/GammaPiPiJpsi.h"])
     .set_constant({ "ECMS" => [:double, 4.660] })
     .note(:multi_energy, "Data taken at multiple energy points: 4.61-4.95 GeV")
     .note(:lepton_separation, "Charged tracks with p<1 GeV/c assigned as pions, p>1 GeV/c as leptons; EMC <0.4 GeV for muons, >1.0 GeV for electrons")
     .note(:jpsi_signal, "M(l+l-) in [3.08,3.12] GeV/c^2; sideband [3.02,3.06] U [3.14,3.18]")
     .note(:isr_psi3686, "Significant psi(3686) peak observed from ISR gamma events")
     .note(:no_x3872_signal, "No obvious X(3872) signal observed; upper limit set on sigma_gamma/sigma_omega < 0.23 at 90% CL")

sel_c = Selection.new

sel_c.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 700
  angle_to_track 10.0
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=1"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :pion, against: [:kaon]
  identify :electron, against: [:pion]
  identify :muon, against: [:pion]
}
.kinematic_fit([:pip, :pim, :ep, :em, :gamma]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_c.with_decay_card(decay_card_gamma_pipi_jpsi).apply(sel_c)
alg_c.execute_on(data_points + incMC_points + [exMC_gamma_pipi_jpsi])