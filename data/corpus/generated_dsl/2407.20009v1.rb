# frozen_string_literal: true

# ============================================================
# Dataset preparation
# ============================================================

# Real data: eight energy points spanning 4.700 – 4.946 GeV
data_4700 = DatasetManager.real_data.find("706_4700")
data_4720 = DatasetManager.real_data.find("706_4720")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")
data_points = [data_4700, data_4720, data_4740, data_4750,
               data_4780, data_4840, data_4914, data_4946]

# Matching inclusive MC samples
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_4720 = DatasetManager.inclusive_mc.find("706_4720")
incMC_4740 = DatasetManager.inclusive_mc.find("707_4740")
incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")
incMC_points = [incMC_4700, incMC_4720, incMC_4740, incMC_4750,
                incMC_4780, incMC_4840, incMC_4914, incMC_4946]

# Signal decay card:
#   e+e- -> K+ K- psi(2S), psi(2S) -> J/psi pi+ pi-, J/psi -> e+ e-
#   RecID map: 0=psi(4260), 1=K+, 2=K-, 3=psi(2S), 4=J/psi, 5=pi+, 6=pi-, 7=e+, 8=e-
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 K+ K- psi(2S) PHSP;
  Enddecay

  Decay psi(2S)
  1.0000 J/psi pi+ pi- PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# 500k-event exclusive signal MC generated at each energy point (common card)
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_KKpsi2S"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMCs.each { |m| m.save_to_config(format: :yaml, file_path: "temp_kksig") }

# All algorithms run on real data + inclusive MC + exclusive signal MC
all_datasets = data_points + incMC_points + exMCs

# ============================================================
# Approach (i): tag K+K- (+ J/psi), signal from recoil mass against K+K-
# ============================================================
alg_i = Algorithm.new("KKpsi2SApprI")
alg_i.set_header(["KKpsi2SApprIAlg/KKpsi2SApprI.h"])
     .set_constant({ "ECMS" => [:double, 4.700] })
     .set_alias({ "std::vector<double>" => "Vdouble" })
alg_i.note(:pid_correction_method,
           "electrons identified by E/p > 0.8 and muons by E/p < 0.4; the J/psi -> mu+mu- "
           "channel additionally requires >=1 muon with hits in >=3 MUC layers")

sel_i = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"
    nChrn     ">=1"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    nkp ">=1"
    nkm ">=1"
  }
  # K+ (recID 1) and K- (recID 2) are tagged; the recoil against K+K- is the psi(2S).
  # Signal region of the psi(2S) recoil mass: 3.67 - 3.71 GeV/c^2.
  .partial_rec([1, 2]) {
    require_recoil_mass 3.67, 3.71
  }
alg_i.with_decay_card(decay_card_signal).apply(sel_i)

# ============================================================
# Approach (ii): tag a single K + psi(2S) -> pi+pi- J/psi (J/psi -> mu+mu-),
#               1C fit to the missing kaon mass (chi2_1C < 50)
# ============================================================
alg_ii = Algorithm.new("KKpsi2SApprII")
alg_ii.set_header(["KKpsi2SApprIIAlg/KKpsi2SApprII.h"])
      .set_constant({ "ECMS" => [:double, 4.700] })
      .set_alias({ "std::vector<double>" => "Vdouble" })
alg_ii.note(:pid_correction_method,
            "electrons identified by E/p > 0.8 and muons by E/p < 0.4; the "
            "J/psi -> mu+mu- selection requires >=1 muon with hits in >=3 MUC layers")

sel_ii = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"
    nChrn     ">=1"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    nkp  "==1"
    npip "==1"
    npim "==1"
    nmup "==1"
    nmum "==1"
  }
  # One K (kp) is tagged; the opposite kaon (km) is left missing and constrained
  # by the 1C missing-kaon fit. J/psi window on mu+mu-: (3.05, 3.15) GeV/c^2.
  .kinematic_fit([:kp, :pip, :pim, :mup, :mum]) {
    nominal
    miss_track_of :km
    invariant_mass_of(:mup, :mum).between(3.05, 3.15)
    chi2_cut 50
  }
alg_ii.with_decay_card(decay_card_signal).apply(sel_ii)

# ============================================================
# Approach (iii): tag K+K- with psi(2S) -> l+l- , 4C fit with chi2 < 200
# ============================================================
alg_iii = Algorithm.new("KKpsi2SApprIII")
alg_iii.set_header(["KKpsi2SApprIIIAlg/KKpsi2SApprIII.h"])
       .set_constant({ "ECMS" => [:double, 4.700] })
       .set_alias({ "std::vector<double>" => "Vdouble" })
alg_iii.note(:pid_correction_method,
             "electrons identified by E/p > 0.8 and muons by E/p < 0.4; the lepton pair is "
             "two opposite-charge tracks with momentum p > 1.0 GeV/c")

sel_iii = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"
    nChrn     ">=1"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    nkp "==1"
    nkm "==1"
    nlp "==1"
    nlm "==1"
  }
  # Fully reconstructed final state K+K- l+l- ; psi(2S) window on l+l-:
  # (3.631, 3.726) GeV/c^2.
  .kinematic_fit([:kp, :km, :lp, :lm]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:lp, :lm).between(3.631, 3.726)
    chi2_cut 200
  }
alg_iii.with_decay_card(decay_card_signal).apply(sel_iii)

# ============================================================
# Approach (iv): tag a single K with psi(2S) -> l+l- (J/psi -> mu+mu-),
#                1C missing-kaon fit (chi2_1C < 15)
# ============================================================
alg_iv = Algorithm.new("KKpsi2SApprIV")
alg_iv.set_header(["KKpsi2SApprIVAlg/KKpsi2SApprIV.h"])
      .set_constant({ "ECMS" => [:double, 4.700] })
      .set_alias({ "std::vector<double>" => "Vdouble" })
alg_iv.note(:pid_correction_method,
            "electrons identified by E/p > 0.8 and muons by E/p < 0.4; the "
            "J/psi -> mu+mu- selection requires >=1 muon with hits in >=3 MUC layers")
      .note(:background_veto,
            "Bhabha-like events vetoed when cos(theta_e+) > 0.85 or cos(theta_e-) < -0.85")

sel_iv = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"
    nChrn     ">=1"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    nkp  "==1"
    nmup "==1"
    nmum "==1"
  }
  # One K (kp) is tagged, the opposite kaon (km) is constrained by the 1C fit;
  # J/psi window on mu+mu-: (3.05, 3.15) GeV/c^2.
  .kinematic_fit([:kp, :mup, :mum]) {
    nominal
    miss_track_of :km
    invariant_mass_of(:mup, :mum).between(3.05, 3.15)
    chi2_cut 15
  }
alg_iv.with_decay_card(decay_card_signal).apply(sel_iv)

# ============================================================
# Execute all four approaches on real data + inclusive MC + signal MC
# ============================================================
root_files_i   = alg_i.execute_on(all_datasets)
root_files_ii  = alg_ii.execute_on(all_datasets)
root_files_iii = alg_iii.execute_on(all_datasets)
root_files_iv  = alg_iv.execute_on(all_datasets)