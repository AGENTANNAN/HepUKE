### Dataset description ###
# 39 real-data energy points spanning 3.51-4.95 GeV (~20 fb^-1) on the BESIII
# energy scan, together with the matching inclusive-MC samples.
# NOTE: the exact 39-point list must be validated against the published list
# (see algorithm.note(:energy_point_list, ...) below).
scan_data  = DatasetManager.real_data.where(cms_energy: { value: 3510.0..4950.0 })
scan_incMC = DatasetManager.inclusive_mc.where(cms_energy: { value: 3510.0..4950.0 })

# Decay card for the signal process e+e- -> p K- K- anti-Xi+ (EvtGen syntax).
# psi(4260) is used as the top mother following the BESIII/KKMC convention for
# non-resonant (continuum-like) final states; the anti-Xi+ is left undetected
# (no decay entry) and is recovered from the missing mass MM(pK-K-).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 p+ K- K- anti-Xi+ PHSP;
    Enddecay

    End
DECAYCARD

# Generate 100k-event exclusive MC at *each* energy point of the scan
# (one ExclusiveMC per related dataset), default cross section.
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_pKKantiXi"   # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "pKKantiXiMM"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})          # psi(4260) top-mother energy
            .set_alias({"std::vector<double>" => "Vdouble"})

# BOSS-side procedures that the DSL cannot express formally
my_algorithm
  .note(:vertex_fit,
        "a common vertex fit of the three charged tracks (p, K-, K-) with chi2 < 100 is
         applied before building MM(pK-K-); the DSL secondary_vertex_fit primitive builds
         a composite particle and cannot express a primary-vertex chi2 selection, so the
         fit is performed in the generated BOSS code")
  .note(:energy_point_list,
        "the 39-point energy list (3.51-4.95 GeV, ~20 fb^-1) must be validated against the
         published list of the paper; the dataset query above selects every real-data sample
         in that energy range and must be cross-checked point by point")
  .note(:ecms_per_point,
        "ECMS is set here to the psi(4260) KKMC top-mother value; in production the beam
         energy is taken per energy point of the scan")

# Build the event selection chain
event_selection = Selection.new
  .select_track {                 # Charged track selection: 1 positive + 2 negative tracks
      cos_theta 0.93              # |cos(theta)| < 0.93
      Vz        10.0              # |Vz| < 10 cm
      Vr        1.0               # Vr < 1 cm
      nChrp     "==1"             # exactly 1 positively charged track
      nChrn     "==2"             # exactly 2 negatively charged tracks
      nNet      "==-1"            # net charge -1 (implied by 1 p + 2 K-)
  }
  .pid(method: :probability) {    # PID by the probability method (default)
      prob_cut 0.001              # PID probability > 0.001
      identify :proton, against: [:kaon, :pion]    # p / p-bar vs K and pi
      identify :kaon,   against: [:pion, :proton]  # K+ / K- vs pi and p
      nprp "==1"                  # demand exactly 1 proton
      nkm  "==2"                  # demand exactly 2 K-
      # Any further identified track beyond the required 1p + 2K- is vetoed by the
      # track-count constraints above (and by nChrp/nChrn in select_track).
  }
  # Partial reconstruction: p (recID 1), K- (recID 2), K- (recID 3) are reconstructed;
  # recID 0 is the top mother (psi(4260)) and recID 4 (anti-Xi+) is left undetected,
  # so the anti-Xi+ is inferred from the recoil / missing mass MM(pK-K-).
  .partial_rec([1, 2, 3]) {
      require_recoil_mass 1.28, 1.38   # MM(pK-K-) in [1.28, 1.38] GeV (Xi mass region)
  }

# Generate the complete algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on the full scan (real data + inclusive MC + signal MC at every point)
root_files = my_algorithm.execute_on(scan_data + scan_incMC + exMCs_signal)