# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) real data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC

# Decay card for psi(2S) -> Sigma+ anti-Sigma- omega
decay_card_omega = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Sigma+ anti-Sigma- omega PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.0000 anti-p- pi0 PHSP;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for psi(2S) -> Sigma+ anti-Sigma- phi
decay_card_phi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Sigma+ anti-Sigma- phi PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay anti-Sigma-
  1.0000 anti-p- pi0 PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K- VSS;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 200k-event exclusive MC for each signal mode
exMC_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_SigmaSigmaBar_omega"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_omega
  config.cross_section   = :default
end

exMC_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_SigmaSigmaBar_phi"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_phi
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ============================ Mode I: psi(2S) -> Sigma+ anti-Sigma- omega ============================
alg_name_omega = "SigmaSigmaBarOmega"
alg_omega = Algorithm.new(alg_name_omega)
alg_omega.set_header(["#{alg_name_omega}Alg/#{alg_name_omega}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})           # 3.686 GeV centre-of-mass energy
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_omega = Selection.new
sel_omega.select_track {            # charged-track selection
    cos_theta 0.93                  # |cos(theta)| < 0.93
    Vz        10.0                  # |Vz| < 10 cm
    Vr        2.0                   # Vr < 2 cm
    nChrp     "==2"                 # two positive tracks
    nChrn     "==2"                 # two negative tracks
    nNet      "==0"                 # net charge zero
  }
  .select_photon {                  # photon selection
    tdc_emc_start     0             # TDC 0
    tdc_emc_end       14            # TDC 14
    angle_to_track    10.0          # at least 10 degrees away from any track
    energyThreshold_b 0.025         # 25 MeV barrel threshold
    energyThreshold_e 0.050         # 50 MeV endcap threshold
    nGam              ">=6"         # >= 6 photons (3 pi0)
  }
  .pid(method: :probability) {      # proton / anti-proton identification
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ and anti-p- vs K and pi
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])     # remove identified (anti-)protons
  .pid(method: :probability) {      # remaining tracks are pions
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  # Nominal 7C fit: 4C energy-momentum + three pi0 mass constraints on the gamma-gamma pairs
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # pi0 (Sigma+)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # pi0 (anti-Sigma-)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # pi0 (omega)
    chi2_cut 45
  }
  # Competing 4C fit of the 5-gamma hypothesis (stores chi2 for wrong-photon-count veto in ROOT)
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }
  # Competing 4C fit of the 7-gamma hypothesis (stores chi2 for wrong-photon-count veto in ROOT)
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }

alg_omega.with_decay_card(decay_card_omega).apply(sel_omega)

# ============================ Mode II: psi(2S) -> Sigma+ anti-Sigma- phi ============================
alg_name_phi = "SigmaSigmaBarPhi"
alg_phi = Algorithm.new(alg_name_phi)
alg_phi.set_header(["#{alg_name_phi}Alg/#{alg_name_phi}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_phi = Selection.new
sel_phi.select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        2.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"       # >= 2 photons
  }
  .pid(method: :probability) {    # identify p, anti-p and exactly one K+, one K-
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion, :proton]
    nprp ">=1"
    nprm ">=1"
    nkp  "==1"
    nkm  "==1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn, :kp <= :chrgp, :km <= :chrgn])
  # 1C Kalman fit: reconstruct one pi0 from a gamma-gamma pair (chi2 < 20)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  }
  # Partial reconstruction: the second pi0 is undetected and inferred from the recoil
  # (recoil = P4_cms - Sigma+ - phi - anti-p- = the missing pi0)
  # recIDs: 1 = Sigma+, 3 = phi, 8 = anti-p-  (9 = the missing pi0)
  .partial_rec([1, 3, 8]) {
    best_combination_by_mass :Sigma+, 1.189   # Sigma+ nominal mass (GeV/c^2)
    best_combination_by_mass :phi,    1.019   # phi nominal mass (GeV/c^2)
    require_recoil_mass 0.10, 0.17            # missing pi0 recoil-mass window (GeV/c^2)
  }

alg_phi.with_decay_card(decay_card_phi).apply(sel_phi)

### Execute on datasets ###
root_files_omega = alg_omega.execute_on([psip_data, psip_incMC, exMC_omega])
root_files_phi   = alg_phi.execute_on([psip_data, psip_incMC, exMC_phi])