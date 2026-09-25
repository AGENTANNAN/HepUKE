# Paper: 2202.00759v2 — Search for X(3872) → π0χc0, π+π−χc0, π0π0χc0
# Data: BOSS 703 energy scan 4.178–4.288 GeV (9.9 fb−1 total)
# Key features: 3 signal channels × 5 χc0 hadronic decay modes, normalization channel
# with lepton PID, create_exclusive_mc_for over energy points
# Rule T1: each independent decay mode → separate Algorithm

### Dataset preparation — BOSS 703 energy scan points ###
data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")
data_4237 = DatasetManager.real_data.find("703_4237")
data_4245 = DatasetManager.real_data.find("703_4245")
data_4246 = DatasetManager.real_data.find("703_4246")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4270 = DatasetManager.real_data.find("703_4270")
data_4280 = DatasetManager.real_data.find("703_4280")

data_points = [data_4180, data_4190, data_4200, data_4210, data_4220,
               data_4230, data_4237, data_4245, data_4246, data_4260,
               data_4270, data_4280]

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4237 = DatasetManager.inclusive_mc.find("703_4237")
incMC_4245 = DatasetManager.inclusive_mc.find("703_4245")
incMC_4246 = DatasetManager.inclusive_mc.find("703_4246")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4270 = DatasetManager.inclusive_mc.find("703_4270")
incMC_4280 = DatasetManager.inclusive_mc.find("703_4280")

incMC_points = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220,
                incMC_4230, incMC_4237, incMC_4245, incMC_4246, incMC_4260,
                incMC_4270, incMC_4280]

### Decay card — common to all channels (KKMC + psi(4260) top mother) ###
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma X_3872 PHSP;
    Enddecay

    Decay X_3872
    1.0000 pi0 chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC — one card over all energy scan points ###
sig_mc = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_X3872"
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ==============================================================================
# Normalization channel: X(3872) → π+π−J/ψ, J/ψ → ℓ+ℓ− (lepton PID)
# ==============================================================================
alg_norm = Algorithm.new("X3872Norm")
alg_norm.set_header(["X3872NormAlg/X3872Norm.h"])
        .set_constant({"ECMS" => [:double, 4.230]})

sel_norm = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
    nNet "==0"
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    tdc_emc_start 0
    tdc_emc_end 14
  }
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.85
    identify :pion, against: [:kaon]
    npip ">=1"
    npim ">=1"
    nlp "==1"
    nlm "==1"
  }
  .remove([:lp <= :chrgp])
  .remove([:lm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pim, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_norm.with_decay_card(decay_card).apply(sel_norm)
root_files_norm = alg_norm.execute_on(data_points + incMC_points + sig_mc)

# ==============================================================================
# Signal channel 1: X(3872) → π0 χc0  — 5 χc0 decay modes
# ==============================================================================

# χc0 → π+π− (2 charged tracks, ≥2 photons for π0)
alg_s1_pp = Algorithm.new("X3872Pi0ChiC0_PiPi")
alg_s1_pp.set_header(["X3872Pi0ChiC0_PiPiAlg/X3872Pi0ChiC0_PiPi.h"])
         .set_constant({"ECMS" => [:double, 4.230]})
sel_s1_pp = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=1"; nChrn ">=1"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=2" }
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }
alg_s1_pp.with_decay_card(decay_card).apply(sel_s1_pp)
root_files_s1_pp = alg_s1_pp.execute_on(data_points + incMC_points + sig_mc)

# χc0 → K+K− (2 kaon tracks, ≥2 photons for π0)
alg_s1_kk = Algorithm.new("X3872Pi0ChiC0_KK")
alg_s1_kk.set_header(["X3872Pi0ChiC0_KKAlg/X3872Pi0ChiC0_KK.h"])
         .set_constant({"ECMS" => [:double, 4.230]})
sel_s1_kk = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=1"; nChrn ">=1"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=2" }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp "==1"; nkm "==1"
  }
  .kinematic_fit([:kp, :km, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }
alg_s1_kk.with_decay_card(decay_card).apply(sel_s1_kk)
root_files_s1_kk = alg_s1_kk.execute_on(data_points + incMC_points + sig_mc)

# χc0 → π+π−π+π− (4 pions, ≥2 photons for π0)
alg_s1_4pi = Algorithm.new("X3872Pi0ChiC0_4pi")
alg_s1_4pi.set_header(["X3872Pi0ChiC0_4piAlg/X3872Pi0ChiC0_4pi.h"])
          .set_constant({"ECMS" => [:double, 4.230]})
sel_s1_4pi = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=2"; nChrn ">=2"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=2" }
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pip, :pim, :pim, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }
alg_s1_4pi.with_decay_card(decay_card).apply(sel_s1_4pi)
root_files_s1_4pi = alg_s1_4pi.execute_on(data_points + incMC_points + sig_mc)

# χc0 → π+π−K+K− (2 pions, 2 kaons, ≥2 photons for π0)
alg_s1_2pi2k = Algorithm.new("X3872Pi0ChiC0_2pi2K")
alg_s1_2pi2k.set_header(["X3872Pi0ChiC0_2pi2KAlg/X3872Pi0ChiC0_2pi2K.h"])
           .set_constant({"ECMS" => [:double, 4.230]})
sel_s1_2pi2k = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=2"; nChrn ">=2"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=2" }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp "==1"; nkm "==1"
  }
  .remove([:kp <= :chrgp])
  .remove([:km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pim, :kp, :km, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }
alg_s1_2pi2k.with_decay_card(decay_card).apply(sel_s1_2pi2k)
root_files_s1_2pi2k = alg_s1_2pi2k.execute_on(data_points + incMC_points + sig_mc)

# χc0 → π+π−π0π0 (2 pions, ≥4 photons, 2π0)
alg_s1_2pi2pi0 = Algorithm.new("X3872Pi0ChiC0_2pi2pi0")
alg_s1_2pi2pi0.set_header(["X3872Pi0ChiC0_2pi2pi0Alg/X3872Pi0ChiC0_2pi2pi0.h"])
             .set_constant({"ECMS" => [:double, 4.230]})
sel_s1_2pi2pi0 = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=1"; nChrn ">=1"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=6" }
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:pip, :pim, :pi0, :pi0, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }
alg_s1_2pi2pi0.with_decay_card(decay_card).apply(sel_s1_2pi2pi0)
root_files_s1_2pi2pi0 = alg_s1_2pi2pi0.execute_on(data_points + incMC_points + sig_mc)

# ==============================================================================
# Signal channel 2: X(3872) → π+π− χc0  — 5 χc0 decay modes
# ==============================================================================

# χc0 → π+π− (4 pions total: signal π+π− + χc0 π+π−)
alg_s2_pp = Algorithm.new("X3872PiPiChiC0_PiPi")
alg_s2_pp.set_header(["X3872PiPiChiC0_PiPiAlg/X3872PiPiChiC0_PiPi.h"])
         .set_constant({"ECMS" => [:double, 4.230]})
sel_s2_pp = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=2"; nChrn ">=2"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14 }
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pim, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_s2_pp.with_decay_card(decay_card).apply(sel_s2_pp)
root_files_s2_pp = alg_s2_pp.execute_on(data_points + incMC_points + sig_mc)

# χc0 → K+K− (2 pions + 2 kaons)
alg_s2_kk = Algorithm.new("X3872PiPiChiC0_KK")
alg_s2_kk.set_header(["X3872PiPiChiC0_KKAlg/X3872PiPiChiC0_KK.h"])
         .set_constant({"ECMS" => [:double, 4.230]})
sel_s2_kk = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=2"; nChrn ">=2"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14 }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp "==1"; nkm "==1"
  }
  .remove([:kp <= :chrgp])
  .remove([:km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pim, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_s2_kk.with_decay_card(decay_card).apply(sel_s2_kk)
root_files_s2_kk = alg_s2_kk.execute_on(data_points + incMC_points + sig_mc)

# χc0 → π+π−π+π− (6 pions total: 2 signal + 4 χc0)
alg_s2_4pi = Algorithm.new("X3872PiPiChiC0_4pi")
alg_s2_4pi.set_header(["X3872PiPiChiC0_4piAlg/X3872PiPiChiC0_4pi.h"])
          .set_constant({"ECMS" => [:double, 4.230]})
sel_s2_4pi = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=3"; nChrn ">=3"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14 }
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pim, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_s2_4pi.with_decay_card(decay_card).apply(sel_s2_4pi)
root_files_s2_4pi = alg_s2_4pi.execute_on(data_points + incMC_points + sig_mc)

# χc0 → π+π−K+K− (4 pions + 2 kaons)
alg_s2_2pi2k = Algorithm.new("X3872PiPiChiC0_2pi2K")
alg_s2_2pi2k.set_header(["X3872PiPiChiC0_2pi2KAlg/X3872PiPiChiC0_2pi2K.h"])
           .set_constant({"ECMS" => [:double, 4.230]})
sel_s2_2pi2k = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=3"; nChrn ">=3"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14 }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp "==1"; nkm "==1"
  }
  .remove([:kp <= :chrgp])
  .remove([:km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pim, :pip, :pim, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_s2_2pi2k.with_decay_card(decay_card).apply(sel_s2_2pi2k)
root_files_s2_2pi2k = alg_s2_2pi2k.execute_on(data_points + incMC_points + sig_mc)

# χc0 → π+π−π0π0 (4 pions + ≥4 photons + 2π0)
alg_s2_2pi2pi0 = Algorithm.new("X3872PiPiChiC0_2pi2pi0")
alg_s2_2pi2pi0.set_header(["X3872PiPiChiC0_2pi2pi0Alg/X3872PiPiChiC0_2pi2pi0.h"])
             .set_constant({"ECMS" => [:double, 4.230]})
sel_s2_2pi2pi0 = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=2"; nChrn ">=2"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=4" }
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:pip, :pim, :pip, :pim, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_s2_2pi2pi0.with_decay_card(decay_card).apply(sel_s2_2pi2pi0)
root_files_s2_2pi2pi0 = alg_s2_2pi2pi0.execute_on(data_points + incMC_points + sig_mc)

# ==============================================================================
# Signal channel 3: X(3872) → π0π0 χc0  — 5 χc0 decay modes
# ==============================================================================

# χc0 → π+π− (2 pions + ≥4 photons, 2π0 from signal + 2γ for π0)
alg_s3_pp = Algorithm.new("X3872Pi0Pi0ChiC0_PiPi")
alg_s3_pp.set_header(["X3872Pi0Pi0ChiC0_PiPiAlg/X3872Pi0Pi0ChiC0_PiPi.h"])
         .set_constant({"ECMS" => [:double, 4.230]})
sel_s3_pp = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=1"; nChrn ">=1"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=4" }
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }
alg_s3_pp.with_decay_card(decay_card).apply(sel_s3_pp)
root_files_s3_pp = alg_s3_pp.execute_on(data_points + incMC_points + sig_mc)

# χc0 → K+K− (2 kaons + ≥4 photons, 2π0 from signal)
alg_s3_kk = Algorithm.new("X3872Pi0Pi0ChiC0_KK")
alg_s3_kk.set_header(["X3872Pi0Pi0ChiC0_KKAlg/X3872Pi0Pi0ChiC0_KK.h"])
         .set_constant({"ECMS" => [:double, 4.230]})
sel_s3_kk = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=1"; nChrn ">=1"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=4" }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp "==1"; nkm "==1"
  }
  .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }
alg_s3_kk.with_decay_card(decay_card).apply(sel_s3_kk)
root_files_s3_kk = alg_s3_kk.execute_on(data_points + incMC_points + sig_mc)

# χc0 → π+π−π+π− (4 pions + ≥4 photons, 2π0 from signal)
alg_s3_4pi = Algorithm.new("X3872Pi0Pi0ChiC0_4pi")
alg_s3_4pi.set_header(["X3872Pi0Pi0ChiC0_4piAlg/X3872Pi0Pi0ChiC0_4pi.h"])
          .set_constant({"ECMS" => [:double, 4.230]})
sel_s3_4pi = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=2"; nChrn ">=2"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=4" }
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pip, :pim, :pim, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }
alg_s3_4pi.with_decay_card(decay_card).apply(sel_s3_4pi)
root_files_s3_4pi = alg_s3_4pi.execute_on(data_points + incMC_points + sig_mc)

# χc0 → π+π−K+K− (2 pions, 2 kaons, ≥4 photons, 2π0 from signal)
alg_s3_2pi2k = Algorithm.new("X3872Pi0Pi0ChiC0_2pi2K")
alg_s3_2pi2k.set_header(["X3872Pi0Pi0ChiC0_2pi2KAlg/X3872Pi0Pi0ChiC0_2pi2K.h"])
           .set_constant({"ECMS" => [:double, 4.230]})
sel_s3_2pi2k = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=2"; nChrn ">=2"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=4" }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp "==1"; nkm "==1"
  }
  .remove([:kp <= :chrgp])
  .remove([:km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:pip, :pim, :kp, :km, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }
alg_s3_2pi2k.with_decay_card(decay_card).apply(sel_s3_2pi2k)
root_files_s3_2pi2k = alg_s3_2pi2k.execute_on(data_points + incMC_points + sig_mc)

# χc0 → π+π−π0π0 (2 pions + ≥8 photons, 2π0 signal + 2π0 χc0)
alg_s3_2pi2pi0 = Algorithm.new("X3872Pi0Pi0ChiC0_2pi2pi0")
alg_s3_2pi2pi0.set_header(["X3872Pi0Pi0ChiC0_2pi2pi0Alg/X3872Pi0Pi0ChiC0_2pi2pi0.h"])
             .set_constant({"ECMS" => [:double, 4.230]})
sel_s3_2pi2pi0 = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp ">=1"; nChrn ">=1"; nNet "==0" }
  .select_photon { energyThreshold_b 0.025; energyThreshold_e 0.050; angle_to_track 10.0;
                    tdc_emc_start 0; tdc_emc_end 14; nGam ">=8" }
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=4"
  }
  .kinematic_fit([:pip, :pim, :pi0, :pi0, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_s3_2pi2pi0.with_decay_card(decay_card).apply(sel_s3_2pi2pi0)
root_files_s3_2pi2pi0 = alg_s3_2pi2pi0.execute_on(data_points + incMC_points + sig_mc)

# Note: χ2/DOF optimization values per channel (Table IV) applied at ROOT stage.
# Simultaneous fit over 5 χc0 modes performed in ROOT analysis.
# Post-fit cuts: χc0 mass window ±25 MeV, X(3872) fit range [3.75, 4.0] GeV/c2.
# J/ψ mass window ±20 MeV/c2 for normalization channel, cosθ_ππ < 0.98,
# η/η' veto cuts for normalization channel applied at ROOT stage.