# frozen_string_literal: true

### Dataset preparation ###
# Full set of real XYZ data points in sqrt(s) ~ 3.81 - 4.95 GeV, together with the
# matching inclusive MC samples (sample name convention: [BOSS]_[CMS energy in MeV]).
xyz_sample_names = %w[
  703_4009 703_4180 703_4190 703_4190scan 703_4200 703_4210 703_4210scan
  703_4220 703_4220scan 703_4230 703_4230scan 703_4237 703_4246 703_4260
  703_4270 703_4280 703_4360 703_4420 703_4600
  705_4130 705_4160
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]

xyz_data  = xyz_sample_names.map { |name| DatasetManager.real_data.find(name) }
xyz_incMC = xyz_sample_names.map { |name| DatasetManager.inclusive_mc.find(name) }

# Mode I decay card: e+e- -> D0 anti-D0, single tag D0 -> K- pi+ pi+ pi- (charge
# conjugate anti-D0 -> K+ pi- pi- pi+ included). The charmed meson opposite the
# tag is left undecayed in the selection (inferred from the recoil mass).
decay_card_D0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- pi- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# Mode II decay card: e+e- -> D+ D-, single tag D+ -> K- pi+ pi+ (charge conjugate
# D- -> K+ pi- pi- included). The recoil D- is not reconstructed.
decay_card_Dp = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K- pi+ pi+ PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# 50k-event exclusive MC for each of the two signal modes, one sample per energy point
exMC_D0 = DatasetManager.create_exclusive_mc_for(xyz_data) do |config|
  config.sample_name   = "exmc_xyz_D0toKPiPiPiPi"
  config.events        = 50_000
  config.decay_card    = decay_card_D0
  config.cross_section = :default
end

exMC_Dp = DatasetManager.create_exclusive_mc_for(xyz_data) do |config|
  config.sample_name   = "exmc_xyz_DplusToKPiPi"
  config.events        = 50_000
  config.decay_card    = decay_card_Dp
  config.cross_section = :default
end

### Event selection (BOSS) ###
# ---------------- Mode I: D0 single tag ----------------
alg_D0 = Algorithm.new("D0STag")
alg_D0.set_header(["D0STagAlg/D0STag.h"])
      .set_constant("ECMS" => [:double, 4.26])   # per-point sqrt(s) taken from each dataset

sel_D0 = Selection.new
sel_D0.select_track {                       # charged track selection
          cos_theta 0.93                    # |cos(theta)| < 0.93
          Vz        10.0                    # |Vz| < 10 cm
          Vr        1.0                     # Vr < 1 cm
          nChrp     ">=2"                   # at least two positive tracks (D0 tag)
          nChrn     ">=2"                   # at least two negative tracks (D0 tag)
        }
        .pid(method: :probability) {        # probability PID
          prob_cut 0.001                    # probability > 0.001
          identify :kaon, against: [:pion, :proton]   # K+ and K- vs pi and p
          identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K and p
          nkm  ">=1"                        # at least one K-
          npip ">=2"                        # at least two pi+
          npim ">=1"                        # at least one pi- (Mode I only)
        }
        # Partial reconstruction: pick the K- pi+ pi+ pi- combination closest to the
        # nominal D0 mass and require the recoil mass to be in the D0 window
        .partial_rec([1]) {                 # recID 1 = D0 (daughters auto-expanded)
          best_combination_by_mass :D0, 1.86484   # nominal D0 mass (GeV/c^2)
          require_recoil_mass 1.80, 1.95          # recoil (anti-D0) mass window
        }

alg_D0.note(:tag_mass_window,
              "the ~14 MeV window around the nominal D0 mass (1.86484 GeV/c^2) cannot " \
              "be enforced by partial_rec, which only selects the combination closest " \
              "to the nominal mass; the hard window is applied downstream")
alg_D0.with_decay_card(decay_card_D0).apply(sel_D0)

# ---------------- Mode II: D+ single tag ----------------
alg_Dp = Algorithm.new("DpSTag")
alg_Dp.set_header(["DpSTagAlg/DpSTag.h"])
      .set_constant("ECMS" => [:double, 4.26])   # per-point sqrt(s) taken from each dataset

sel_Dp = Selection.new
sel_Dp.select_track {                       # charged track selection
          cos_theta 0.93                    # |cos(theta)| < 0.93
          Vz        10.0                    # |Vz| < 10 cm
          Vr        1.0                     # Vr < 1 cm
          nChrp     ">=2"                   # at least two positive tracks (D+ tag)
          nChrn     ">=1"                   # at least one negative track (D+ tag)
        }
        .pid(method: :probability) {        # probability PID
          prob_cut 0.001                    # probability > 0.001
          identify :kaon, against: [:pion, :proton]   # K+ and K- vs pi and p
          identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K and p
          nkm  ">=1"                        # at least one K-
          npip ">=2"                        # at least two pi+
        }
        # Partial reconstruction: pick the K- pi+ pi+ combination closest to the
        # nominal D+ mass and require the recoil mass to be in the D0 window
        .partial_rec([1]) {                 # recID 1 = D+ (daughters auto-expanded)
          best_combination_by_mass :Dp, 1.86966   # nominal D+ mass (GeV/c^2)
          require_recoil_mass 1.80, 1.95          # recoil (D-) mass window
        }

alg_Dp.note(:tag_mass_window,
              "the ~16 MeV window around the nominal D+ mass (1.86966 GeV/c^2) cannot " \
              "be enforced by partial_rec, which only selects the combination closest " \
              "to the nominal mass; the hard window is applied downstream")
alg_Dp.with_decay_card(decay_card_Dp).apply(sel_Dp)

# ---------------- Execute on real data, inclusive MC and signal MC ----------------
root_files_modeI  = alg_D0.execute_on(xyz_data + xyz_incMC + exMC_D0)
root_files_modeII = alg_Dp.execute_on(xyz_data + xyz_incMC + exMC_Dp)