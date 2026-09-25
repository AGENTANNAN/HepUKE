# frozen_string_literal: true

### Dataset preparation ###
# e+e- -> p pbar for the proton electromagnetic form factors.
# 22 R-scan c.m. energy points from 2.00 to 3.08 GeV (BOSS 713, round08).
rscan_points = [2000, 2050, 2100, 2125, 2150, 2175, 2200, 2232, 2309, 2386,
                2396, 2500, 2644, 2646, 2700, 2800, 2900, 2950, 2981, 3000,
                3020, 3080]

# Real data and the corresponding inclusive MC sample at each scan point.
rscan_data  = rscan_points.map { |e| DatasetManager.real_data.find("713_Rscan_#{e}") }
rscan_incMC = rscan_points.map { |e| DatasetManager.inclusive_mc.find("713_Rscan_#{e}") }

# ConExc decay card for e+e- -> p pbar.  The generator models ISR (up to second
# order) and the vacuum polarisation by itself.  `Particle vpho` is intentionally
# omitted: the DSL injects the per-point c.m. energy for every scan point.
decay_card_ppbar = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc 8;
  Enddecay
  End
DECAYCARD

# Signal-exclusive MC: 500k events per c.m. energy point (one sample per point,
# the c.m. energy being injected for each scan point).
exMC_ppbar = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "exmc_ppbar_rscan"
  config.events        = 500_000
  config.decay_card    = decay_card_ppbar
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name  = "PpbarFormFactor"
ppbar_alg = Algorithm.new(alg_name)
ppbar_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({
           "ECMS"    => [:double, 3.08],  # nominal c.m. energy; the per-scan-point value is used at run time
           "P_MEAN"  => [:double, 2.0],   # per-scan-point fitted mean momentum (GeV/c)
           "P_SIGMA" => [:double, 0.05]   # per-scan-point fitted momentum width (GeV/c)
         })
         .set_alias({"std::vector<double>" => "Vdouble"})

# BOSS-side procedures that have no DSL counterpart are preserved as notes.
ppbar_alg.note(:opening_angle_cut, "opening angle between p and pbar in the c.m. frame required to be > 170 deg at 2.00/2.05 GeV, > 175 deg for 2.100-2.309 GeV and > 178 deg for 2.386-3.080 GeV")
         .note(:cosmic_rejection, "cosmic-ray events rejected by requiring |T1 - T2| < 4 ns, where T1 and T2 are the TOF times of the two charged tracks")
         .note(:background_veto, "E/p cut applied above 2.150 GeV to suppress the Bhabha (e+e- -> e+e-) background")
         .note(:vertex_chi2_cut, "the two charged tracks are required to form a common vertex with chi2 < 100")

event_selection = Selection.new
event_selection.select_track {
                 cos_theta 0.93    # |cos(theta)| < 0.93 for each track
                 Vz        100.0   # |Vz| < 100 cm
                 Vr        10.0    # Vr < 10 cm
                 nChrp     "==1"   # exactly one positive track
                 nChrn     "==1"   # exactly one negative track
                 nNet      "==0"   # net charge zero
               }
               .pid(method: :probability) {
                 prob_cut 0.001    # PID probability > 0.001
                 identify :proton, against: [:kaon, :pion]  # p and pbar vs K and pi
                 nprp "==1"        # one proton
                 nprm "==1"        # one antiproton
               }
               .remove([:prp <= :chrgp, :prm <= :chrgn])   # remove the identified p/pbar from the track lists
               # Asymmetric momentum window (p_mean - 4*sigma) < p < (p_mean + 3*sigma),
               # fitted independently at each c.m. energy.
               .remove(:prp) { condition "three_momentum_of(:prp) < (P_MEAN - 4.0*P_SIGMA) || three_momentum_of(:prp) > (P_MEAN + 3.0*P_SIGMA)" }
               .remove(:prm) { condition "three_momentum_of(:prm) < (P_MEAN - 4.0*P_SIGMA) || three_momentum_of(:prm) > (P_MEAN + 3.0*P_SIGMA)" }
               .kinematic_fit([:prp, :prm]) {
                 nominal
                 vertex_fit([0, 1])         # p and pbar must come from a common vertex
                 constrain_four_momentum    # 4C constraint of the ppbar system to the c.m. energy
                 chi2_cut 200
               }

ppbar_alg.with_decay_card(decay_card_ppbar).apply(event_selection)

# Real data, inclusive MC and the ConExc signal MC all go through the same chain.
root_files = ppbar_alg.execute_on(rscan_data + rscan_incMC + exMC_ppbar)