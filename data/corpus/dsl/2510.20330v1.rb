# Precision Measurement of Ds*+ - Ds+ Mass Difference
# BESIII Collaboration, arXiv:2510.20330v1
# Data: 3.19 fb^-1 at sqrt(s)=4.178 GeV (BOSS 703_4180)
# Signal: Ds*+ -> Ds+ pi0, Ds+ -> K+ K- pi+
# Calibration: D*+ -> D+ pi0, D+ -> K- pi+ pi+

### Dataset preparation ###
data_4180 = DatasetManager.real_data.find("703_4180")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

# Signal decay card: Ds*+ -> Ds+ pi0, Ds+ -> K+ K- pi+
decay_card_Ds = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s-         PHSP;
    Enddecay

    Decay D_s*+
    1.0000 D_s+ pi0           VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.0000 K+ K- pi+          PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi-          PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma        PHSP;
    Enddecay

    End
DECAYCARD

exMC_Ds = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "DsStar_Ds_mass_diff_signal"
  config.related_dataset = data_4180
  config.events = 100000
  config.decay_card = decay_card_Ds
  config.cross_section = :default
end

# Calibration decay card: D*+ -> D+ pi0, D+ -> K- pi+ pi+
decay_card_D = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D*+ D-             PHSP;
    Enddecay

    Decay D*+
    1.0000 D+ pi0             VSP_PWAVE;
    Enddecay

    Decay D+
    1.0000 K- pi+ pi+         PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-         PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma        PHSP;
    Enddecay

    End
DECAYCARD

exMC_D = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "DStar_D_mass_diff_calibration"
  config.related_dataset = data_4180
  config.events = 100000
  config.decay_card = decay_card_D
  config.cross_section = :default
end

### Event selection: Signal channel (Ds*+ -> Ds+ pi0, Ds+ -> K+ K- pi+) ###
alg_Ds = Algorithm.new("DsStarMassDiff")
alg_Ds.set_header(["DsStarMassDiffAlg/DsStarMassDiff.h"])
       .set_constant({"ECMS" => [:double, 4.178]})

sel_Ds = Selection.new
sel_Ds.select_track {
         cos_theta 0.93
         Vz 10.0
         Vr 1.0
         nChrp ">=2"
         nChrn ">=1"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :kaon, against: [:pion]
         identify :pion, against: [:kaon]
         nkp ">=1"
         nkm ">=1"
         npip ">=1"
       }
       .remove([:kp <= :chrgp])
       .remove([:km <= :chrgn])
       .remove([:pip <= :chrgp])
       .select_photon {
         tdc_emc_start 0
         tdc_emc_end 14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         nGam ">=2"
       }
       .kalman_kinematic_fit([:gamma, :gamma]) do
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
         chi2_cut 25
         npi0 ">=1"
       end
       .kinematic_fit([:kp, :km, :pip, :pi0]) do
         nominal
         constrain_four_momentum
         chi2_cut 200
       end
# Post-kinematic-fit steps (ROOT):
# - Ds+ mass window: |M(K+K-pi+) - m_Ds+| < 12 MeV/c^2
# - p(pi0) < 100 MeV/c for soft pi0 from Ds*+ decay
# - Recoil mass cuts against tag Ds-
# - Delta_M = M(K+K-pi+pi0) - M(K+K-pi+) mass difference fit
# - CB + Bifurcated Gaussian + Gaussian composite PDF fit

alg_Ds.with_decay_card(decay_card_Ds).apply(sel_Ds)


### Event selection: Calibration channel (D*+ -> D+ pi0, D+ -> K- pi+ pi+) ###
alg_D = Algorithm.new("DStarMassDiffCalib")
alg_D.set_header(["DStarMassDiffCalibAlg/DStarMassDiffCalib.h"])
      .set_constant({"ECMS" => [:double, 4.178]})

sel_D = Selection.new
sel_D.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=2"
        nChrn ">=1"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion]
        identify :pion, against: [:kaon]
        nkm ">=1"
        npip ">=2"
      }
      .remove([:km <= :chrgn])
      .remove([:pip <= :chrgp])
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=2"
      }
      .kalman_kinematic_fit([:gamma, :gamma]) do
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      end
      .kinematic_fit([:km, :pip, :pip, :pi0]) do
        nominal
        constrain_four_momentum
        chi2_cut 200
      end
# Post-kinematic-fit steps (ROOT):
# - D+ mass window: |M(K-pi+pi+) - m_D+| < 12 MeV/c^2
# - p(pi0) < 100 MeV/c for soft pi0 from D*+ decay
# - Recoil mass cuts
# - Delta_M = M(K-pi+pi+pi0) - M(K-pi+pi+) mass difference
# - Calibration of Ds*+ mass difference using well-known D*+ mass difference

alg_D.with_decay_card(decay_card_D).apply(sel_D)

### Execute ###
root_files_Ds = alg_Ds.execute_on([data_4180, incMC_4180, exMC_Ds])
root_files_D = alg_D.execute_on([data_4180, incMC_4180, exMC_D])