### Dataset description ###
# Real data and inclusive MC at the two energy points (4.914 and 4.946 GeV, BOSS 7.0.7)
data_4914  = DatasetManager.real_data.find("707_4914")
data_4946  = DatasetManager.real_data.find("707_4946")
incmc_4914 = DatasetManager.inclusive_mc.find("707_4914")
incmc_4946 = DatasetManager.inclusive_mc.find("707_4946")

# Decay card: e+e- -> phi chi_c1(3872) -> K+K- pi+pi- l+l-   (electron channel)
# Top mother set to psi(4260) following the BESIII KKMC convention for a direct
# continuum-like production of the signal final state.
decay_card_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0000 phi chi_c1(3872) PHSP;
  Enddecay

  Decay chi_c1(3872)
  1.0000 rho0 J/psi PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K- VSS;
  Enddecay

  Decay rho0
  1.0000 pi+ pi- VSS;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Decay card: same chain, muon channel J/psi -> mu+mu-
decay_card_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 phi chi_c1(3872) PHSP;
  Enddecay

  Decay chi_c1(3872)
  1.0000 rho0 J/psi PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K- VSS;
  Enddecay

  Decay rho0
  1.0000 pi+ pi- VSS;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# 50k-event exclusive MC for the full signal chain, generated at both energy points
exMCs_ee = DatasetManager.create_exclusive_mc_for([data_4914, data_4946]) do |config|
  config.sample_name   = "phi_chic1_3872_ee"
  config.events        = 50000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

exMCs_mumu = DatasetManager.create_exclusive_mc_for([data_4914, data_4946]) do |config|
  config.sample_name   = "phi_chic1_3872_mumu"
  config.events        = 50000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ================= Class 1: 6-track  K+K-pi+pi-l+l-  =================
alg_6trk = Algorithm.new("PhiChiC1SixTrack")
alg_6trk.set_header(["PhiChiC1SixTrackAlg/PhiChiC1SixTrack.h"])
        .set_constant({"ECMS" => [:double, 4.914]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_6trk = Selection.new
sel_6trk.select_track {                 # charged-track selection
          cos_theta 0.93                # |cos(theta)| < 0.93
          Vz        10.0                # |Vz| < 10 cm
          Vr        1.0                 # Vr < 1 cm
          nChrp     "==3"               # three positive tracks
          nChrn     "==3"               # three negative tracks
          nNet      "==0"               # zero net charge
        }
        .select_photon {                # photon selection (no multiplicity cut)
          tdc_emc_start     0           # TDC 0-14
          tdc_emc_end       14
          angle_to_track    10.0        # >= 10 deg from any charged track
          energyThreshold_b 0.025       # EMC barrel energy > 25 MeV
          energyThreshold_e 0.050       # EMC endcap energy > 50 MeV
        }
        .pid(method: :probability) {    # probability PID
          prob_cut 0.001                # probability cut for K/pi separation
          # tracks with p > 1.0 GeV are leptons; e if EMC energy > 0.6 GeV, else mu
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                         treat_as_electron_if_energy_above: 0.6
          identify :kaon, against: [:pion]
          identify :pion, against: [:kaon]
        }
        # 4C kinematic fit of K+K-pi+pi-l+l-
        .kinematic_fit([:kp, :km, :pip, :pim, :lp, :lm]) {
          nominal
          constrain_four_momentum
          chi2_cut 150
        }

alg_6trk.with_decay_card(decay_card_ee).apply(sel_6trk)
alg_6trk.execute_on([data_4914, data_4946, incmc_4914, incmc_4946] + exMCs_ee + exMCs_mumu)

# ================= Class 2: 5-track (one missing kaon) =================
alg_5trk = Algorithm.new("PhiChiC1FiveTrack")
alg_5trk.set_header(["PhiChiC1FiveTrackAlg/PhiChiC1FiveTrack.h"])
        .set_constant({"ECMS" => [:double, 4.914]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_5trk = Selection.new
sel_5trk.select_track {                 # charged-track selection
          cos_theta 0.93                # |cos(theta)| < 0.93
          Vz        10.0                # |Vz| < 10 cm
          Vr        1.0                 # Vr < 1 cm
          nChrp     ">=2"               # at least two positive tracks
          nChrn     ">=2"               # at least two negative tracks
          nTot      "==5"               # five tracks in total
          nNet      "==0"               # zero net charge
        }
        .select_photon {                # photon selection (no multiplicity cut)
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.050
        }
        .pid(method: :probability) {    # probability PID
          prob_cut 0.001
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                         treat_as_electron_if_energy_above: 0.6
          identify :kaon, against: [:pion]
          identify :pion, against: [:kaon]
        }
        # 1C fit, hypothesis: K- is the missing kaon (seen kaon is K+)
        .kinematic_fit([:kp, :pip, :pim, :lp, :lm]) {
          nominal
          miss_track_of :km
          constrain_four_momentum
          chi2_cut 20
        }
        # 1C fit, alternative hypothesis: K+ is the missing kaon (seen kaon is K-);
        # the two stored chi2 values allow the better hypothesis to be selected.
        .kinematic_fit([:km, :pip, :pim, :lp, :lm]) {
          miss_track_of :kp
          constrain_four_momentum
          chi2_cut 20
        }

alg_5trk.with_decay_card(decay_card_ee).apply(sel_5trk)
alg_5trk.execute_on([data_4914, data_4946, incmc_4914, incmc_4946] + exMCs_ee + exMCs_mumu)