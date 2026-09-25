# frozen_string_literal: true
# e+e- -> Xi- anti-Xi+ (single-baryon tag), sqrt(s) = 4.009 - 4.600 GeV
# BOSS part: dataset preparation + event selection up to the partial reconstruction.

### Dataset preparation ###
# 15 XYZ (703-1) energy points spanning 4.009 - 4.600 GeV, total ~11 fb^-1
scan_points = %w[
  703_4009 703_4180 703_4190 703_4200 703_4210
  703_4220 703_4230 703_4237 703_4246 703_4260
  703_4270 703_4280 703_4360 703_4420 703_4600
]
# Corresponding real data and inclusive MC sample at each energy point
real_data    = scan_points.map { |name| DatasetManager.real_data.find(name) }
inclusive_mc = scan_points.map { |name| DatasetManager.inclusive_mc.find(name) }

# Decay card for the signal process (EvtGen format).
# e+e- -> Xi- anti-Xi+ with Xi- -> pi- Lambda0, Lambda0 -> p+ pi- and the charge conjugates.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Xi- anti-Xi+ PHSP;
    Enddecay

    Decay Xi-
    1.0000 pi- Lambda0 PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 pi+ anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC generated for every energy point (same card / cross section).
exMC_signal = DatasetManager.create_exclusive_mc_for(real_data) do |config|
  config.sample_name   = "exmc_xi_minus_xi_bar_plus"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "XiMinusTag"
xi_algorithm = Algorithm.new(alg_name)
xi_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 4.009] }) # representative value; per-point energy comes from the dataset
            .set_alias({ "std::vector<double>" => "Vdouble" })

# BOSS-side procedures that the DSL cannot express are kept as notes.
xi_algorithm
  .note(:track_vertex_cut,
        "no Vz / Vr (vertex) cut is applied in the charged-track selection; only |cos(theta)| < 0.93
         and the charged-track multiplicities are required")
  .note(:lambda_selection_cuts,
        "Lambda -> p pi- secondary-vertex fit requires chi2 < 500 for 3 d.o.f.,
         |M(p pi-) - 1.115683 GeV/c^2| < 5 MeV/c^2 and decay length > 0; the vertex chi2 / mass window /
         flight-length cuts are not expressible in the current DSL")
  .note(:xi_minus_selection_cuts,
        "Xi- built from pi- Lambda by mass-difference minimisation, selecting the combination closest to the
         nominal Xi- mass 1.32171 GeV/c^2; |M(pi- Lambda) - 1.32171 GeV/c^2| < 10 MeV/c^2 and decay length > 0
         are required; these windows / flight-length cuts are not expressible in the current DSL")

event_selection = Selection.new
event_selection
  .select_track {          # charged-track selection (no Vz/Vr cut)
    cos_theta 0.93         # |cos(theta)| < 0.93
    nChrp     ">=1"        # at least one positively charged track
    nChrn     ">=2"        # at least two negatively charged tracks
  }
  .pid(method: :probability) {   # particle identification: mutual p / K / pi separation
    prob_cut 0.001               # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]     # p+ and p- (charge-conjugation shorthand)
    identify :kaon,   against: [:pion, :proton]   # K+ and K-
    identify :pion,   against: [:kaon, :proton]   # pi+ and pi-
    nprp ">=1"                    # at least one proton
    npim ">=2"                    # at least two pi-
  }
  # remove the identified K+- and p/pbar from the charged-track lists
  .remove([:kp <= :chrgp, :km <= :chrgn, :prp <= :chrgp, :prm <= :chrgn])
  # the remaining positive / negative charged tracks are taken as pi+ / pi-
  .assign({ :chrgp => :pip, :chrgn => :pim })
  # Lambda -> p pi- via a secondary-vertex fit (best combination by mass difference)
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Xi- -> pi- Lambda via a secondary-vertex fit (best combination by mass difference)
  .secondary_vertex_fit([:pim, :Lambda]) {
    build_virtual_particle(:Xi_minus).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Partial reconstruction: anti-Xi+ (recID 2, daughters expanded automatically) is left missing;
  # the recoil mass against the tagged Xi- is required to be in [1.2, 1.5] GeV/c^2.
  # No global kinematic fit is performed (partial reconstruction replaces it).
  .partial_miss([2]) {
    require_recoil_mass 1.2, 1.5
  }

# Attach the decay card and render the full selection into the BOSS algorithm.
xi_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Run on real data, inclusive MC and the per-point exclusive signal MC.
root_files = xi_algorithm.execute_on(real_data + inclusive_mc + exMC_signal)