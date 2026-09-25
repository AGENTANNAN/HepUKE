# =============================================================================
# BOSS part : dataset preparation + event selection
#   e+e- -> Ds*+ DsJ- ,  Ds*+ -> gamma Ds+ ,  Ds+ -> K+ K- pi+ ,
#   DsJ- -> pi+ pi- pi0    for DsJ- = Ds0*(2317)- , Ds1(2460)- , Ds1(2536)-
#   data : 4.600 - 4.700 GeV, seven energy points
# =============================================================================

### ------------------------------- Datasets ------------------------------- ###
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.600 GeV
data_4610 = DatasetManager.real_data.find("706_4610")   # 4.612 GeV
data_4620 = DatasetManager.real_data.find("706_4620")   # 4.628 GeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4.641 GeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4.661 GeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4.682 GeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4.699 GeV
data_points = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]

# Matching inclusive MC samples for the same seven energy points
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_samples = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

### ------------------------------ Decay cards ------------------------------ ###
# DsJ- = Ds0*(2317)-
decay_card_Ds0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s0*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s+
  1.000 K+ K- pi+ PHSP;
  Enddecay

  Decay D_s0*-
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# DsJ- = Ds1(2460)-
decay_card_Ds1 = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s1- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s+
  1.000 K+ K- pi+ PHSP;
  Enddecay

  Decay D_s1-
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# DsJ- = Ds1(2536)-
decay_card_Ds1p = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s1prime- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s+
  1.000 K+ K- pi+ PHSP;
  Enddecay

  Decay D_s1prime-
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### ---------------------------- Exclusive MC ---------------------------- ###
# 500k events per energy point, for each of the three DsJ- signal hypotheses
exMC_2317 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_DsStarDsJ2317"   # auto-suffixed per energy point
  config.events        = 500_000
  config.decay_card    = decay_card_Ds0
  config.cross_section = :default
end

exMC_2460 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_DsStarDsJ2460"
  config.events        = 500_000
  config.decay_card    = decay_card_Ds1
  config.cross_section = :default
end

exMC_2536 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_DsStarDsJ2536"
  config.events        = 500_000
  config.decay_card    = decay_card_Ds1p
  config.cross_section = :default
end

### ------------------------------ Algorithms ------------------------------ ###
# ECMS is only a nominal constant here : the per-point beam energy of the scan
# is supplied by each dataset / job.  One algorithm per DsJ- decay card.
dsplus_submode_note =
  "Ds+ candidates are accepted only if AT LEAST ONE of the two resonant submodes " \
  "is satisfied : Ds+ -> phi pi+  with |M(K+K-) - m(phi)| < 9 MeV/c2,  or " \
  "Ds+ -> anti-K*0 K+  with |M(K-pi+) - m(K*0)| < 84 MeV/c2.  The DSL expresses a " \
  "single mass window (invariant_mass_of(...).within / .out_of) but not the logical OR " \
  "of two windows, so this disjunction is applied afterwards in the ROOT analysis."

alg_2317 = Algorithm.new("DsStarDsJ2317")
alg_2317.set_header(["DsStarDsJ2317Alg/DsStarDsJ2317.h"])
        .set_constant({"ECMS" => [:double, 4.65]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:dsplus_submode_selection, dsplus_submode_note)

alg_2460 = Algorithm.new("DsStarDsJ2460")
alg_2460.set_header(["DsStarDsJ2460Alg/DsStarDsJ2460.h"])
        .set_constant({"ECMS" => [:double, 4.65]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:dsplus_submode_selection, dsplus_submode_note)

alg_2536 = Algorithm.new("DsStarDsJ2536")
alg_2536.set_header(["DsStarDsJ2536Alg/DsStarDsJ2536.h"])
        .set_constant({"ECMS" => [:double, 4.65]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:dsplus_submode_selection, dsplus_submode_note)

### ---------------------------- Event selection ---------------------------- ###
# The three DsJ- hypotheses share the same reconstructed final state
# (gamma K+ K- pi+, with the DsJ- undetected) and the same selection chain.
event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93        # |cos(theta)| < 0.93
    Vz        10.0        # |Vz| < 10 cm
    Vr        1.0         # Vr < 1 cm
    nChrp     "==2"       # two positive tracks (K+ and pi+)
    nChrn     "==1"       # one negative track (K-)
    nNet      "==1"       # net charge +1
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025   # 25 MeV in the EMC barrel
    energyThreshold_e 0.050   # 50 MeV in the EMC endcap
    angle_to_track    20.0    # > 20 degrees away from any charged track
    nGam              ">=1"   # at least one photon (from Ds*+ -> gamma Ds+)
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]   # K+ and K- (separated from pions)
    identify :pion, against: [:kaon]   # pi+ and pi- (separated from kaons)
    nkp  "==1"                         # exactly one K+
    nkm  "==1"                         # exactly one K-
    npip "==1"                         # exactly one pi+
  }
  # 2C kinematic fit on gamma K+ K- pi+ :
  #   m(K+ K- pi+)          = m(Ds+)  = 1.96834 GeV/c2
  #   m(gamma K+ K- pi+)    = m(Ds*+) = 2.1122  GeV/c2
  # No 4C energy-momentum constraint is applied (2 constraints only).
  .kinematic_fit([:gamma, :kp, :km, :pip]) {
    nominal
    invariant_mass_of(:kp, :km, :pip).constrain_to_nominal_mass_of(:"D_s+")
    invariant_mass_of(:gamma, :kp, :km, :pip).constrain_to_nominal_mass_of(:"D_s*+")
    chi2_cut 10
  }
  # Partial reconstruction : Ds*+ is built from gamma K+ K- pi+ (with its
  # corrected four-momenta); the DsJ- is NOT reconstructed and is inferred from
  # the recoil four-momentum.
  # recIDs of the decay card : 0 psi(4260) | 1 D_s*+ | 2 D_sJ- | 3 gamma |
  #   4 D_s+ | 5 K+ | 6 K- | 7 pi+(Ds+) | 8 pi+(DsJ-) | 9 pi-(DsJ-) |
  #   10 pi0 | 11 gamma | 12 gamma
  .partial_rec([1, 3, 4, 5, 6, 7]) {
    best_combination_by_mass :"D_s*+", 2.1122   # combination closest to nominal m(Ds*+)
    require_recoil_mass 2.0, 2.8                # M(recoil) = M(DsJ-) in 2.0 - 2.8 GeV/c2
  }

### --------------------------- Apply and execute --------------------------- ###
alg_2317.with_decay_card(decay_card_Ds0).apply(event_selection)
alg_2460.with_decay_card(decay_card_Ds1).apply(event_selection.dup)
alg_2536.with_decay_card(decay_card_Ds1p).apply(event_selection.dup)

root_files_2317 = alg_2317.execute_on(data_points + incMC_samples + exMC_2317)
root_files_2460 = alg_2460.execute_on(data_points + incMC_samples + exMC_2460)
root_files_2536 = alg_2536.execute_on(data_points + incMC_samples + exMC_2536)