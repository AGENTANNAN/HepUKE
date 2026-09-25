# ======================================================================
# e+e- -> p pbar on the psi(2S) scan  (continuum, ConExc)
# BOSS part: dataset preparation + full event selection up to the
#            nominal kinematic fit.
# ======================================================================

### ---------------------------------------------------------------
### Dataset preparation
### ---------------------------------------------------------------
# Nine real-data CM energy points of the psi(2S) scan (BOSS 704),
# sqrt(s) = 3.581 - 3.710 GeV, total luminosity ~ 495 pb^-1.
scan_points = [
  "704_3581",   # 3581.5 MeV,  84.6 pb^-1
  "704_3670",   # 3670.2 MeV,  83.6 pb^-1
  "704_3680",   # 3680.1 MeV,  83.1 pb^-1
  "704_3683",   # 3682.8 MeV,  28.2 pb^-1
  "704_3684",   # 3684.2 MeV,  27.8 pb^-1
  "704_3685",   # 3685.3 MeV,  25.3 pb^-1
  "704_3686",   # 3686.5 MeV,  24.5 pb^-1
  "704_3691",   # 3691.4 MeV,  68.6 pb^-1
  "704_3710"    # 3709.8 MeV,  69.3 pb^-1
].map { |sample_name| DatasetManager.real_data.find(sample_name) }

# Signal decay card - ConExc generator (continuum / ISR Born cross section).
# The signal is produced through a single virtual photon with ISR modelled up to
# O(alpha^2) and the measured Born cross section; the p pbar pair itself is
# generated in phase space.  There is therefore NO KKMC / psi(4260) top mother,
# and `Particle vpho <ECMS> 0.0` is deliberately NOT written here: for a
# multi-point energy scan the DSL injects the correct sqrt(s) at every point.
decay_card_ppbar = <<~DECAYCARD
    Decay vpho
    1.0000 ConExc 74110;      # ConExc mode for vpho -> p pbar (phase space, ISR-corrected)
    Enddecay

    End
DECAYCARD

# Matching exclusive MC, one sample per scan point (100k events per point);
# create_exclusive_mc_for re-uses the card/cross-section and injects the correct
# sqrt(s) (Particle vpho) at each energy point.
exMC_ppbar = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_ppbar_conexc"   # auto-suffixed per scan point
  config.events        = 100_000               # 100k events at every point
  config.decay_card    = decay_card_ppbar
  config.cross_section = :default
end

### ---------------------------------------------------------------
### Event selection (BOSS)
### ---------------------------------------------------------------
alg_name = "PpbarContinuum"
ppbar_alg = Algorithm.new(alg_name)
ppbar_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         # Nominal psi(2S) energy of the scan (see note: the per-point sqrt(s)
         # spans 3.581 - 3.710 GeV and is used point-by-point in the fit).
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

# Full selection chain - 1:1 with the described steps.
event_selection = Selection.new
event_selection
  # Charged tracks: barrel MDC only
  .select_track {
      cos_theta 0.80        # |cos(theta)| < 0.80 (barrel)
      Vz        10.0        # |Vz| < 10 cm
      Vr        1.0         # Vr < 1 cm in the transverse plane
      nChrp     "==1"       # exactly one positive track
      nChrn     "==1"       # exactly one negative track
      nNet      "==0"       # net charge zero
  }
  # Particle identification: p / pbar against K and pi
  .pid(method: :probability) {
      prob_cut 0.001                                  # PID probability > 0.001
      identify :proton, against: [:kaon, :pion]       # p and pbar at once
      nprp "==1"                                      # exactly one proton
      nprm "==1"                                      # exactly one anti-proton
  }
  # E/p < 0.5 (raw EMC energy over MDC momentum) for both candidates:
  # removes the e+e- (Bhabha) and mu+mu- backgrounds.  Candidates failing the
  # cut are erased from the proton / anti-proton candidate lists.
  .for_each(:prp) do
      where { eraw_over_p >= 0.5 }   # active subset: candidates failing E/p < 0.5
      remove                         # erase them
  end
  .for_each(:prm) do
      where { eraw_over_p >= 0.5 }
      remove
  end
  # Nominal 4C kinematic fit on the two proton candidates; the corrected
  # four-momenta are the ones stored for the ROOT-side momentum window.
  .kinematic_fit([:prp, :prm]) {
      nominal                  # nominal fit - corrected four-momenta are saved
      constrain_four_momentum  # constrain total four-momentum to the CMS energy
      chi2_cut 200             # loose cut in BOSS; tight cut optimised in ROOT
  }

# BOSS-side procedures that have no dedicated DSL primitive.
ppbar_alg
  .note(:cosmic_veto, "cosmic-ray events are rejected by requiring the two TOF
    flight times to satisfy |T(p) - T(pbar)| < 3 ns; the TOF time difference of
    the p/pbar pair has no dedicated DSL primitive and is applied to the PID'd
    candidates before the kinematic fit")
  .note(:back_to_back_cut, "the opening angle between the p and pbar momenta,
    evaluated in the e+e- CM frame, is required to lie within [178, 180] deg;
    the event-by-event CM boost of the two raw MDC momenta is not expressible
    with the available DSL primitives (applied before the kinematic fit)")
  .note(:per_point_cms_energy, "this is an energy scan: the ECMS constant must
    take the measured sqrt(s) of each scan point (3.5815, 3.6702, 3.6801,
    3.6828, 3.6842, 3.6853, 3.6865, 3.6914, 3.7098 GeV) when running that point,
    since it enters the 4C constraint; the value 3.686 in the spec is only the
    nominal psi(2S) label of the scan")

# Attach the decay card, render the selection, run on data + matching signal MC.
ppbar_alg.with_decay_card(decay_card_ppbar).apply(event_selection)

root_files = ppbar_alg.execute_on(scan_points + exMC_ppbar)