# ================================================================
# e+e- -> p K- Lambda (charge conjugate implied) at six BESIII points:
# 4.008, 4.178, 4.226, 4.258, 4.416, 4.682 GeV
# Aim: spin-parity determination of X(2085)
# ================================================================

### Dataset preparation ###
data_points = [
  DatasetManager.real_data.find("703_4009"),  # 4.008 GeV
  DatasetManager.real_data.find("703_4180"),  # 4.178 GeV
  DatasetManager.real_data.find("703_4230"),  # 4.226 GeV
  DatasetManager.real_data.find("703_4260"),  # 4.258 GeV
  DatasetManager.real_data.find("703_4420"),  # 4.416 GeV
  DatasetManager.real_data.find("706_4680")   # 4.682 GeV
]

incMC_points = [
  DatasetManager.inclusive_mc.find("703_4009"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("706_4680")
]

# Signal decay card: psi(4260) -> p+ K- Lambda (KKMC generator), Lambda -> p pi- (phase space)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 p+ K- Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC per energy point (same card, six datasets)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_pKLambda_X2085"   # auto-suffixed per energy point
  config.events = 100_000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "X2085pKLambda"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.258]})  # CMS energy (set per energy point at run time)
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                       # Charged-track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        20.0                      # |Vz| < 20 cm
    nChrp     ">=2"                     # at least 2 positive tracks (prompt p, Lambda -> p)
    nChrn     ">=2"                     # at least 2 negative tracks (K-, Lambda -> pi-)
    nNet      "==0"                     # net charge zero
  }
  .pid(method: :probability) {          # PID, probability method
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ / pbar
    identify :kaon,   against: [:pion, :proton] # K+ / K-
    identify :pion,   against: [:kaon, :proton] # pi+ / pi-  (Lambda -> pi-)
    nprp ">=2"                          # at least two protons
    nkm  ">=1"                          # at least one K-
  }
  .secondary_vertex_fit([:prp, :pim]) { # Reconstruct Lambda -> p pi- via secondary vertex
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list  # remaining tracks are the prompt p and K-
  }
  .kinematic_fit([:prp, :km, :Lambda]) {        # 4C kinematic fit to p K- Lambda
    nominal
    constrain_four_momentum
    chi2_cut 100
  }

# BOSS-side procedures that have no dedicated DSL method
alg.note(:lambda_vertex_quality, "Lambda -> p pi- secondary vertex required chi2(vertex) < 100 and a decay length > 2 sigma of the vertex resolution")
   .note(:lambda_mass_window, "Lambda candidate required |M(p pi-) - M(Lambda)| < 6 MeV/c^2")
   .note(:prompt_track_quality, "the two tracks remaining after building the Lambda (prompt p and K-) required tighter |dz| < 10 cm and |dr| < 1 cm")
   .note(:background_veto, "events with |cos theta_K| > 0.83 rejected to suppress Bhabha background")
   .note(:multi_energy_ecms, "the identical selection is applied at all six energy points (4.008, 4.178, 4.226, 4.258, 4.416, 4.682 GeV); the CMS-energy constant ECMS must be set per dataset, and the value shown is representative")

alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and exclusive signal MC at all six points
root_files = alg.execute_on(data_points + incMC_points + exMCs_signal)