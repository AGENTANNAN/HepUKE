# Paper 2405.07741v2: Search for chi_c1(3872) -> gamma psi_2(3823) at 4.178-4.278 GeV
# Signal: e+e- -> gamma chi_c1(3872), chi_c1(3872) -> gamma psi_2(3823),
#         psi_2(3823) -> gamma chi_c1, chi_c1 -> gamma J/psi, J/psi -> l+ l-
# Final state: l+ l- + 4 photons

# --- Datasets (BOSS 703, 4.178-4.278 GeV) ---
data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")
data_4237 = DatasetManager.real_data.find("703_4237")
data_4246 = DatasetManager.real_data.find("703_4246")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4270 = DatasetManager.real_data.find("703_4270")
data_4280 = DatasetManager.real_data.find("703_4280")

data_points = [data_4180, data_4190, data_4200, data_4210, data_4220,
               data_4230, data_4237, data_4246, data_4260, data_4270, data_4280]

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4237 = DatasetManager.inclusive_mc.find("703_4237")
incMC_4246 = DatasetManager.inclusive_mc.find("703_4246")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4270 = DatasetManager.inclusive_mc.find("703_4270")
incMC_4280 = DatasetManager.inclusive_mc.find("703_4280")

incMC_points = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220,
                incMC_4230, incMC_4237, incMC_4246, incMC_4260, incMC_4270, incMC_4280]

# --- Decay card: e+e- -> gamma chi_c1(3872) -> four gamma cascade + J/psi -> l+ l- ---
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 gamma chi_c1(3872) PHSP;
    Enddecay

    Decay chi_c1(3872)
    1.000 gamma psi_2(3823) PHSP;
    Enddecay

    Decay psi_2(3823)
    1.000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC ---
exMC = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_X3872_gammaPsi23823"
  config.events = 200000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ==========================================================================
# Algorithm: Search for chi_c1(3872) -> gamma psi_2(3823)
# 2 charged lepton tracks + at least 4 photons
# EMC-based e/mu separation, 4C + 7C kinematic fits
# ==========================================================================
alg = Algorithm.new("X3872ToGammaPsi23823")
alg.set_header(["X3872ToGammaPsi23823Alg/X3872ToGammaPsi23823.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })
   .note(:multi_energy, "Data taken at 11 energy points: 4.178-4.278 GeV; run separately per dataset")
   .note(:lepton_id, "Exactly 2 charged tracks. EMC deposit >= 0.8 GeV for electrons, < 0.4 GeV for muons")
   .note(:jpsi_window, "|M(ll) - m(J/psi)| < 30 MeV/c^2 for J/psi selection")
   .note(:pi0_veto, "All photon pair invariant masses must be > 15 MeV/c^2 from nominal pi0 mass to suppress pi0 pi0 J/psi background")
   .note(:four_c_fit, "4C kinematic fit: constrain 4-momentum of lepton pair + 4 photons to beam; best combination by minimum chi2")
   .note(:seven_c_fit, "7C kinematic fit: 4C + mass constraints on M(ll)=m_J/psi, M(gamma ll)=m_chi_c1, M(gamma gamma ll)=m_psi_2(3823); chi2_7C < 100")
   .note(:photon_assignment, "7C fit resolves which photon comes from each radiative cascade; best combination by minimum chi2_7C")
   .note(:e1_transition, "e+e- -> gamma X(3872) production simulated as E1 transition; chi_c1(3872)->gamma psi_2(3823) and psi_2(3823)->gamma chi_c1 use phase space")
   .note(:result, "No signal observed. Upper limit R < 0.075 at 90% CL on B(chi_c1(3872)->gamma psi_2(3823))/B(chi_c1(3872)->pi+pi-J/psi)")

event_selection = Selection.new

event_selection.select_track {
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
  energyThreshold_e 0.025
  nGam ">=4"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :electron, against: [:pion]
  identify :muon, against: [:pion]
}
.kinematic_fit([:ep, :em, :gamma, :gamma, :gamma, :gamma]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg.with_decay_card(decay_card).apply(event_selection)
alg.execute_on(data_points + incMC_points + [exMC])