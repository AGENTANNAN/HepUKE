# =====================================================================
# BESIII detector-performance measurement at sqrt(s) = 4.61 - 4.95 GeV
# Three independent measurements:
#   (1) CMS energy   : e+e- -> Lambda_c+ Lambda_c-  (partial reconstruction)
#   (2) luminosity   : e+e- -> (gamma) e+e-        (large-angle Bhabha)
#   (3) cross-check  : e+e- -> (gamma) gamma gamma (di-photon)
# No kinematic fit is performed at the BOSS level; the energy and
# luminosity extractions are done offline in ROOT.
# =====================================================================

### Dataset description ###
# --- twelve real-data energy points (BOSS 706: 4.610-4.700 GeV; BOSS 707: 4.740-4.946 GeV) ---
data_points = [
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4914"),
  DatasetManager.real_data.find("707_4946")
]

# --- corresponding inclusive MC samples ---
incMC_points = [
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4914"),
  DatasetManager.inclusive_mc.find("707_4946")
]

# --- decay card for (1): psi(4260) -> Lambda_c+ Lambda_c- , phase-space ---
decay_card_lambdac = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 p+ K- pi+ PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.000 anti-p- K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- decay card for (2): e+e- -> (gamma) e+e- (Bhabha) ---
decay_card_bhabha = <<~DECAYCARD
  Decay psi(4260)
  1.000 e+ e- PHOTOS PHSP;
  Enddecay

  End
DECAYCARD

# --- decay card for (3): e+e- -> (gamma) gamma gamma ---
decay_card_diphoton = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma gamma PHOTOS PHSP;
  Enddecay

  End
DECAYCARD

# Signal MC for the Lambda_c+Lambda_c- channel: one sample per energy point
exMC_lambdac = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_psi4260_lambdacpair"
  config.events        = 100000
  config.decay_card    = decay_card_lambdac
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ---------------------------------------------------------------------
# (1) CMS energy from Lambda_c+ Lambda_c- partial reconstruction
#     Lambda_c+ -> p K- pi+ is reconstructed; the opposite Lambda_c- is
#     inferred from the recoil four-momentum. No kinematic fit.
# ---------------------------------------------------------------------
alg_name_lc = "LambdaCPairEnergy"
alg_lc = Algorithm.new(alg_name_lc)
alg_lc.set_header(["#{alg_name_lc}Alg/#{alg_name_lc}.h"])
      .set_constant({ "ECMS" => [:double, 4.610] })   # beam energy varies over the 12-point scan
      .set_alias({ "std::vector<double>" => "Vdouble" })

sel_lc = Selection.new
sel_lc.select_track {            # common track-quality selection
        cos_theta 0.93           # |cos(theta)| < 0.93
        Vz 10.0                  # |Vz| < 10 cm
        Vr 1.0                   # Vr < 1 cm
      }
      .pid(method: :probability) {   # common PID: separate K, pi and p from one another
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :kaon,   against: [:pion, :proton]
        identify :pion,   against: [:kaon, :proton]
      }
      .partial_rec([1]) {            # reconstruct Lambda_c+ (recID 1) and its daughters; infer the rest from recoil
        best_combination_by_mass :"Lambda_c+", 2.28646   # pick the p K- pi+ combination closest to the Lambda_c mass
      }

alg_lc.note(:ecms_per_point, "The beam energy ECMS varies over the twelve scan points
  (4.610-4.946 GeV); the CMS four-momentum used in the recoil calculation is taken per run
  from the measured beam energy, and the final energy calibration is performed in ROOT.")
alg_lc.with_decay_card(decay_card_lambdac).apply(sel_lc)

# ---------------------------------------------------------------------
# (2) Luminosity from large-angle Bhabha scattering e+e- -> (gamma) e+e-
# ---------------------------------------------------------------------
alg_name_bh = "BhabhaLuminosity"
alg_bh = Algorithm.new(alg_name_bh)
alg_bh.set_header(["#{alg_name_bh}Alg/#{alg_name_bh}.h"])
      .set_constant({ "ECMS" => [:double, 4.610] })
      .set_alias({ "std::vector<double>" => "Vdouble" })

sel_bh = Selection.new
sel_bh.select_track {            # exactly two oppositely charged tracks, EMC barrel acceptance
        cos_theta 0.80           # |cos(theta)| < 0.8 (EMC barrel)
        Vz 10.0
        Vr 1.0
        nChrp "==1"              # exactly one positive track
        nChrn "==1"              # exactly one negative track
        nNet  "==0"              # net charge zero
      }
      .remove(:chrgp) { condition "three_momentum_of(:chrgp) < 2.0" }   # momentum > 2 GeV/c
      .remove(:chrgn) { condition "three_momentum_of(:chrgn) < 2.0" }
      .pid(method: :probability) {     # common PID plus lepton identification
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6
        identify :proton, against: [:kaon, :pion]
        identify :kaon,   against: [:pion, :proton]
        identify :pion,   against: [:kaon, :proton]
        nlp "==1"
        nlm "==1"
      }

alg_bh.note(:emc_saturation_correction, "For high-energy Bhabha electrons the EMC energy
  deposit is used to correct the shower saturation effect before the energy/momentum
  comparison; the correction curve is applied offline in ROOT.")
alg_bh.with_decay_card(decay_card_bhabha).apply(sel_bh)

# ---------------------------------------------------------------------
# (3) Di-photon cross-check e+e- -> (gamma) gamma gamma
# ---------------------------------------------------------------------
alg_name_gg = "DiPhotonCrossCheck"
alg_gg = Algorithm.new(alg_name_gg)
alg_gg.set_header(["#{alg_name_gg}Alg/#{alg_name_gg}.h"])
      .set_constant({ "ECMS" => [:double, 4.610] })
      .set_alias({ "std::vector<double>" => "Vdouble" })

sel_gg = Selection.new
sel_gg.select_track {            # common track-quality selection
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
      }
      .pid(method: :probability) {     # common PID: separate K, pi and p from one another
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :kaon,   against: [:pion, :proton]
        identify :pion,   against: [:kaon, :proton]
      }
      .select_photon {               # two photon candidates in the EMC
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0          # > 10 deg from the nearest charged track
        energyThreshold_b 0.025      # barrel energy threshold
        energyThreshold_e 0.050      # endcap energy threshold
        nGam "==2"                   # exactly two photons
      }

alg_gg.with_decay_card(decay_card_diphoton).apply(sel_gg)

### Execute ###
all_datasets = data_points + incMC_points
root_files_lc = alg_lc.execute_on(all_datasets + exMC_lambdac)
root_files_bh = alg_bh.execute_on(all_datasets)
root_files_gg = alg_gg.execute_on(all_datasets)