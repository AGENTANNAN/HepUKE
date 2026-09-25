# frozen_string_literal: true

### ------------------------------------------------------------------ ###
### Dataset preparation                                               ###
### ------------------------------------------------------------------ ###
# e+e- -> p pbar gamma_ISR studied at seven c.m. energies (total ~7.5 fb^-1):
#   3.773, 4.180, 4.230, 4.260, 4.360, 4.420, 4.600 GeV
data_points = [
  DatasetManager.real_data.find("712_3773"),   # 3.773 GeV
  DatasetManager.real_data.find("703_4180"),   # 4.180 GeV
  DatasetManager.real_data.find("703_4230"),   # 4.230 GeV
  DatasetManager.real_data.find("703_4260"),   # 4.260 GeV
  DatasetManager.real_data.find("703_4360"),   # 4.360 GeV
  DatasetManager.real_data.find("703_4420"),   # 4.420 GeV
  DatasetManager.real_data.find("703_4600"),   # 4.600 GeV
]

# Matching inclusive MC sample at each energy point
incMC_points = [
  DatasetManager.inclusive_mc.find("712_3773"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4600"),
]

# Decay card for the signal channel, KKMC based (psi(4260) top-mother convention),
# psi(4260) -> p+ anti-p- gamma in phase space
decay_card_ppbar_gamma = <<~DECAYCARD
    Decay psi(4260)
    1.000 p+ anti-p- gamma PHSP;
    Enddecay
    End
DECAYCARD

# 200k-event exclusive MC per energy point: identical decay card / cross section /
# statistics, only the related real dataset differs -> one ExclusiveMC per point
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ppbar_gamma_isr"
  config.events        = 200_000
  config.decay_card    = decay_card_ppbar_gamma
  config.cross_section = :default
end

### ------------------------------------------------------------------ ###
### Event selection (BOSS)                                            ###
### ------------------------------------------------------------------ ###
alg_name = "PPbarGammaISR"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.260]})   # compile-time default; per-energy value set in jobOptions
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  # Charged track selection
  .select_track {
    cos_theta 0.93    # |cos(theta)| < 0.93
    Vz        10.0    # |Vz| < 10 cm
    Vr        1.0     # Vr < 1 cm
    nChrp     "==1"   # exactly one positive track
    nChrn     "==1"   # exactly one negative track
    nNet      "==0"   # net charge zero
  }
  # Photon selection
  .select_photon {
    tdc_emc_start     0      # EMC TDC start
    tdc_emc_end       14     # EMC TDC end
    angle_to_track    10.0   # > 10 deg from nearest charged track
    energyThreshold_b 0.025  # E > 25 MeV (barrel)
    energyThreshold_e 0.050  # E > 50 MeV (endcap)
    nGam              ">=1"  # at least one photon
  }
  # Particle identification: one proton and one anti-proton
  .pid(method: :probability) {
    prob_cut   0.001
    identify :proton, against: [:kaon, :pion]  # p+ and anti-p-, separated from K and pi
    nprp       "==1"                            # exactly one proton
    nprm       "==1"                            # exactly one anti-proton
  }
  # E/p < 0.5 for both (anti)proton candidates -> suppress Bhabha (e+e- -> e+e-)
  .for_each(:prp) {
    where { e_over_p > 0.5 }
    remove
  }
  .for_each(:prm) {
    where { e_over_p > 0.5 }
    remove
  }
  # ISR photon must be isolated by > 20 deg from BOTH the proton and anti-proton tracks
  .select_isolated_photon {
    angle_to_prp_track 20.0
    angle_to_prm_track 20.0
    nGam               ">=1"
  }
  # Remove photons below 0.4 GeV ...
  .for_each(:gamma) {
    where { energy < 0.4 }
    remove
  }
  # ... and tag the highest-energy surviving photon as the ISR photon
  .for_each(:gamma) {
    best { maximize { energy } }
    remove
  }
  # 4C kinematic fit on p pbar gamma (no invariant-mass window applied)
  .kinematic_fit([:prp, :prm, :gamma]) {
    nominal
    constrain_four_momentum  # constrain total four-momentum to the CMS energy
    chi2_cut 50
  }

my_algorithm
  .note(:multi_energy_ecms, "the 4C fit CMS energy must match each of the seven scan points; the ECMS constant here is only a compile-time default and is overridden per energy point in the jobOptions")
  .with_decay_card(decay_card_ppbar_gamma)
  .apply(event_selection)

# Run on real data, inclusive MC and exclusive MC at every energy point
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs)