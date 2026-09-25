# Paper 2004.07701v2: Sigma+ spin polarization in J/psi / psi(3686) -> Sigma+ anti-Sigma-
# Two independent datasets: J/psi and psi(2S) are separate analyses
# Sigma+ -> p+ pi0, anti-Sigma- -> anti-p- pi0, pi0 -> gamma gamma
# Proton/anti-proton PID, pi0 via 1C Kalman fit, 4C kinematic fit

decay_card = <<~DECAYCARD
Decay J/psi
1.0 Sigma+ anti-Sigma- PHSP;
Enddecay
Decay psi(2S)
1.0 Sigma+ anti-Sigma- PHSP;
Enddecay
Decay Sigma+
1.0 p+ pi0 PHSP;
Enddecay
Decay anti-Sigma-
1.0 anti-p- pi0 PHSP;
Enddecay
Decay pi0
1.0 gamma gamma PHSP;
Enddecay
End
DECAYCARD

# --- J/psi ---
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

exMC_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_Sigma_SigmaBar_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg_jpsi = Algorithm.new("Jpsi_to_Sigma_SigmaBar")
alg_jpsi.set_header(["Jpsi_to_Sigma_SigmaBarAlg/Jpsi_to_Sigma_SigmaBar.h"])
alg_jpsi.set_constant({ "ECMS" => [:double, 3.097] })

sel_jpsi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 2.0
    nChrp "==1"
    nChrn "==1"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  }
  # Remove protons from charged track lists, assign remaining charged tracks
  .remove([:prp <= :chrgp])
  .remove([:prm <= :chrgn])
  # Reconstruct two pi0 candidates via Kalman fits
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # 4C kinematic fit
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {
    constrain_four_momentum
    chi2_cut 100
    nominal
  }

# The DSL cannot express: best combination by minimizing
# sqrt((M(p pi0) - m_Sigma+)^2 + (M(anti-p pi0) - m_anti-Sigma-)^2)
alg_jpsi.note(:best_combination,
  "Best combination selected by minimizing sqrt((M(p pi0) - m_Sigma+)^2 + (M(anti-p pi0) - m_anti-Sigma-)^2)."
)

alg_jpsi.with_decay_card(decay_card).apply(sel_jpsi)
alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi])

# --- psi(2S) ---
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

exMC_psip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_Sigma_SigmaBar_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg_psip = Algorithm.new("psip_to_Sigma_SigmaBar")
alg_psip.set_header(["psip_to_Sigma_SigmaBarAlg/psip_to_Sigma_SigmaBar.h"])
alg_psip.set_constant({ "ECMS" => [:double, 3.686] })

sel_psip = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 2.0
    nChrp "==1"
    nChrn "==1"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  }
  .remove([:prp <= :chrgp])
  .remove([:prm <= :chrgn])
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {
    constrain_four_momentum
    chi2_cut 100
    nominal
  }

alg_psip.note(:best_combination,
  "Best combination selected by minimizing sqrt((M(p pi0) - m_Sigma+)^2 + (M(anti-p pi0) - m_anti-Sigma-)^2). " \
  "For psi(2S) additional veto: |M(p anti-p) - 3.1 GeV| > 0.05 GeV to suppress psi(2S)->pi0 pi0 J/psi."
)

alg_psip.with_decay_card(decay_card).apply(sel_psip)
alg_psip.execute_on([psip_data, psip_incMC, exMC_psip])