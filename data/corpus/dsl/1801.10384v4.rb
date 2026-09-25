# BOSS Ruby DSL for BESIII paper 1801.10384v4
# Search for Z_s in e+e- → φππ and cross section measurement at √s = 2.125 GeV
# Two independent decay channels:
#   Mode I:  e+e- → φπ+π-  (φ → K+K-)
#   Mode II: e+e- → φπ0π0 (φ → K+K-, π0 → γγ)

### Dataset description ###
data_2125 = DatasetManager.real_data.find("713_Rscan_2125")
incMC_2125 = DatasetManager.inclusive_mc.find("713_Rscan_2125")

# Decay card: e+e- → φπ+π- (φ → K+K-)
decay_card_charged = <<~DECAYCARD
    Decay psi(4260)
    1.000  phi  pi+  pi-  PHSP;
    Enddecay

    Decay phi
    1.000  K+  K-  VSS;
    Enddecay

    End
DECAYCARD

# Decay card: e+e- → φπ0π0 (φ → K+K-, π0 → γγ)
decay_card_neutral = <<~DECAYCARD
    Decay psi(4260)
    1.000  phi  pi0  pi0  PHSP;
    Enddecay

    Decay phi
    1.000  K+  K-  VSS;
    Enddecay

    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for charged channel
exMC_charged = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_2125_phipipi"
  config.related_dataset = data_2125
  config.events = 100000
  config.decay_card = decay_card_charged
  config.cross_section = :default
end

# Exclusive MC for neutral channel
exMC_neutral = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_2125_phipi0pi0"
  config.related_dataset = data_2125
  config.events = 100000
  config.decay_card = decay_card_neutral
  config.cross_section = :default
end

### Event selection - Mode I: e+e- → φπ+π- (φ → K+K-) ###
# Selection uses a 1C kinematic fit with one kaon missing to maximise efficiency.
# Charged tracks: |cosθ| < 0.93, |Vz| < 10 cm, Vr < 1 cm, N_tot ≥ 3.
# PID: ≥1 kaon, 2 pions of opposite charge.
# 1C fit: detected K+π+π- with missing K-, χ² < 10.
# Paper also handles the complementary case (K- detected, K+ missing).

alg_modeI = Algorithm.new("PhiPiPiCharged")
alg_modeI.set_header(["PhiPiPiChargedAlg/PhiPiPiCharged.h"])
          .set_constant({"ECMS" => [:double, 2.125]})

sel_modeI = Selection.new
sel_modeI.select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        ">=3"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"
    npip ">=1"; npim ">=1"
  end
  .kinematic_fit([:kp, :pip, :pim]) do
    nominal
    miss_track_of :km
    constrain_four_momentum
    chi2_cut 10
  end
  .kinematic_fit([:km, :pip, :pim]) do
    miss_track_of :kp
    constrain_four_momentum
  end

alg_modeI
  .note(:helix_correction, "helix parameter correction applied to charged tracks before kinematic fit; efficiency difference between with/without estimated by re-running BOSS selection")
  .note(:pid_correction_method, "PID efficiency difference between data and MC: 1% per pion, 3% per kaon; MDC tracking efficiency difference < 1.5% per charged track")
  .with_decay_card(decay_card_charged)
  .apply(sel_modeI)

### Event selection - Mode II: e+e- → φπ0π0 (φ → K+K-, π0 → γγ) ###
# Charged tracks: |cosθ| < 0.93, |Vz| < 10 cm, Vr < 1 cm, N_tot ≥ 1.
# Photons: ≥ 4, E > 25 MeV (barrel) / 50 MeV (endcap), angle to track > 10°, 0 ≤ t_emc ≤ 700 ns.
# PID: ≥ 1 kaon.
# π0 reconstruction: 2C Kalman fit (γγ → π0 mass constraint), χ² < 25, ≥ 2 π0.
# 1C kinematic fit: detected K+π0π0 with missing K-, χ² < 20.

alg_modeII = Algorithm.new("PhiPiPiNeutral")
alg_modeII.set_header(["PhiPiPiNeutralAlg/PhiPiPiNeutral.h"])
           .set_constant({"ECMS" => [:double, 2.125]})

sel_modeII = Selection.new
sel_modeII.select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        ">=1"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=4"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"
  end
  .kalman_kinematic_fit([:gamma, :gamma, :gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  end
  .kinematic_fit([:kp, :pi0, :pi0]) do
    nominal
    miss_track_of :km
    constrain_four_momentum
    chi2_cut 20
  end
  .kinematic_fit([:km, :pi0, :pi0]) do
    miss_track_of :kp
    constrain_four_momentum
  end

alg_modeII
  .note(:helix_correction, "helix parameter correction applied to charged tracks before kinematic fit")
  .note(:pid_correction_method, "PID efficiency difference between data and MC: 1% per pion, 3% per kaon; photon selection efficiency difference < 1% per photon (4% total for 4 photons)")
  .note(:efficiency_curve, "π0 mass window requirement: M(γγ) within ±20 MeV/c² of nominal π0 mass, applied at ROOT level")
  .with_decay_card(decay_card_neutral)
  .apply(sel_modeII)

### Execute ###
alg_modeI.execute_on([data_2125, incMC_2125, exMC_charged])
alg_modeII.execute_on([data_2125, incMC_2125, exMC_neutral])