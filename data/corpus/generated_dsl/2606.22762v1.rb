# =============================================================================
# Born cross section for e+e- -> p pbar
#   47 c.m. energies, 3.510 - 4.946 GeV, 26 fb^-1 of real data
#   matching inclusive MC; continuum signal MC generated with the ConExc
#   (ISR) generator at every scan point.
# =============================================================================

### Dataset description ###
# --- 47 real-data scan points (3.510 - 4.946 GeV) ---------------------------
names_703 = %w[chi_c1_scan_4 chi_c1_scan_5 3810 3872 3900 4009 4090 4180 4190
               4200 4210 4220 4230 4237 4245 4246 4260 4270 4280 4310 4360 4390
               4420 4470 4530 4575 4600]
names_705 = %w[4130 4160 4290 4315 4340 4380 4400 4440]
names_706 = %w[4610 4620 4640 4660 4680 4700]
names_707 = %w[4740 4750 4780 4840 4914 4946]

scan_points = names_703.map { |n| DatasetManager.real_data.find("703_#{n}") } +
              names_705.map { |n| DatasetManager.real_data.find("705_#{n}") } +
              names_706.map { |n| DatasetManager.real_data.find("706_#{n}") } +
              names_707.map { |n| DatasetManager.real_data.find("707_#{n}") }

# --- matching inclusive MC over the same energy range ------------------------
scan_inclusive_mc = DatasetManager.inclusive_mc.where(cms_energy: { value: 3510..4946 })

# --- ConExc decay card for the continuum process e+e- -> p pbar ---------------
# The DSL detects the literal ConExc token and switches to the no-KKMC template;
# for a multi-energy scan it injects `Particle vpho <ECMS> 0.0` for each scan
# point, so the Particle vpho line is deliberately omitted here.
continuum_decay_card_ppbar = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 8;
    Enddecay
    End
DECAYCARD

# --- 100k-event exclusive MC at each of the 47 scan points -------------------
# (identical card / cross section; the correct sqrt(s) is injected per point)
exMC_ppbar = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_ppbar_continuum"
  config.events        = 100_000
  config.decay_card    = continuum_decay_card_ppbar
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PpbarScan"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
# No single ECMS constant is declared: the c.m. energy differs at each of the
# 47 scan points and is injected per point by the ConExc scan mechanism.

event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.80      # |cos(theta)| < 0.80
    Vz        10.0      # |Vz| < 10 cm
    Vr        1.0       # Vr < 1 cm
    nChrp     "==1"     # exactly one positive charged track
    nChrn     "==1"     # exactly one negative charged track
    nNet      "==0"     # net charge zero
  }
  .pid(method: :probability) {
    prob_cut 0.001                              # PID likelihood > 0.001
    identify :proton, against: [:kaon, :pion]   # proton hypothesis dominant over K and pi (p+ and pbar)
    nprp "==1"                                  # exactly one proton
    nprm "==1"                                  # exactly one antiproton
  }
  .kinematic_fit([:prp, :prm]) {
    nominal                 # nominal fit -> corrected four-momenta are stored
    constrain_four_momentum # 4C energy-momentum constraint
    chi2_cut 200            # chi2 < 200
  }
  # The remaining background suppression (E/p(p+) < 0.5, p-pbar opening angle
  # > 3.1 rad, MUC hit depth < 40 cm, |p(pbar) - p_expected| < 3 sigma_p) acts on
  # kinematic-fit-corrected quantities and is therefore applied at the ROOT stage.

my_algorithm
  .note(:isr_correction, "ISR radiative corrections obtained by iterative weighting of the ConExc p pbar signal MC separately at each c.m. energy point")
  .note(:efficiency_curve, "detection efficiency determined per c.m. energy by iterative MC weighting; not a flat constant")

my_algorithm.with_decay_card(continuum_decay_card_ppbar).apply(event_selection)
root_files = my_algorithm.execute_on(scan_points + scan_inclusive_mc + exMC_ppbar)