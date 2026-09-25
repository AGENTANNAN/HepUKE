# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Real data and inclusive MC at the two energy points
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # J/psi inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) (3.686 GeV) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # psi(3686) inclusive MC

# Decay card: e+e- -> Sigma+ Sigma- ; Sigma+ -> p pi0 ; Sigma- -> anti-p pi0 ; pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Sigma+ anti-Sigma- PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0 PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: one sample per energy point, same decay chain, 1e6 events each
exMCs_signal = DatasetManager.create_exclusive_mc_for([jpsi_data, psip_data]) do |config|
  config.sample_name   = "sigma_sigmabar_signal"
  config.events        = 1_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMCs_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
# Identical selection chain used at both energy points
common_selection = Selection.new
  .select_track {
      cos_theta 0.93        # |cos(theta)| < 0.93
      Vz        10.0        # |Vz| < 10 cm
      Vr        2.0         # Vr < 2 cm
      nChrp     "==1"       # exactly one positively charged track
      nChrn     "==1"       # exactly one negatively charged track
      nNet      "==0"       # net charge zero
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025    # 25 MeV in the barrel
      energyThreshold_e 0.050    # 50 MeV in the endcap
      angle_to_track    10.0     # >= 10 degrees from any charged track
      nGam              ">=2"    # at least two photons
  }
  .pid(method: :probability) {
      prob_cut 0.001                                 # PID probability > 0.001
      identify :proton, against: [:kaon, :pion]     # p+ and anti-p- separated from K and pi
      nprp "==1"                                    # exactly one proton
      nprm "==1"                                    # exactly one anti-proton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])         # remove identified protons from generic charged lists
  .select_isolated_photon {
      angle_to_prp_track 20.0    # >= 20 degrees from the prompt proton track
      angle_to_prm_track 20.0    # >= 20 degrees from the prompt anti-proton track
      nGam ">=2"                 # at least two isolated photons
  }
  .assign({:chrgp => :pip, :chrgn => :pim})         # remaining tracks assumed to be pi+/pi-
  .kalman_kinematic_fit([:gamma, :gamma]) {         # reconstruct pi0 from gamma gamma
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200
      npi0 ">=1"                                    # at least one pi0 candidate
  }
  .kinematic_fit([:prp, :pi0, :prm, :pi0]) {        # 2C kinematic fit to p pi0 anti-p pi0
      nominal
      constrain_four_momentum                       # four-momentum conservation
      chi2_cut 30
  }

### Algorithm at the J/psi energy (3.097 GeV) ###
alg_jpsi = Algorithm.new("SigmaSigmaBarJpsi")
alg_jpsi.set_header(["SigmaSigmaBarJpsiAlg/SigmaSigmaBarJpsi.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        .set_alias({"std::vector<double>" => "Vdouble"})
alg_jpsi.note(:dca_cut, "proton/anti-proton distance of closest approach to the interaction point (DCA > 0.34 cm) is noted as a potential background-suppression criterion but is NOT applied in this selection.")
alg_jpsi.with_decay_card(decay_card_signal).apply(common_selection.dup)
root_files_jpsi = alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMCs_signal[0]])

### Algorithm at the psi(3686) energy (3.686 GeV) ###
alg_psip = Algorithm.new("SigmaSigmaBarPsip")
alg_psip.set_header(["SigmaSigmaBarPsipAlg/SigmaSigmaBarPsip.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})
alg_psip.note(:dca_cut, "proton/anti-proton distance of closest approach to the interaction point (DCA > 0.34 cm) is noted as a potential background-suppression criterion but is NOT applied in this selection.")
alg_psip.with_decay_card(decay_card_signal).apply(common_selection.dup)
root_files_psip = alg_psip.execute_on([psip_data, psip_incMC, exMCs_signal[1]])