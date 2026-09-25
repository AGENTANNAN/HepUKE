# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# ---- Real data at the seven centre-of-mass energies ----
data_3773 = DatasetManager.real_data.find("712_3773")   # sqrt(s) = 3.773 GeV
data_4008 = DatasetManager.real_data.find("703_4009")   # sqrt(s) ~ 4.008 GeV
data_4226 = DatasetManager.real_data.find("703_4230")   # sqrt(s) ~ 4.226 GeV
data_4258 = DatasetManager.real_data.find("703_4260")   # sqrt(s) ~ 4.258 GeV
data_4358 = DatasetManager.real_data.find("703_4360")   # sqrt(s) ~ 4.358 GeV
data_4416 = DatasetManager.real_data.find("703_4420")   # sqrt(s) ~ 4.416 GeV
data_4600 = DatasetManager.real_data.find("703_4600")   # sqrt(s) = 4.600 GeV

# ---- Inclusive MC at the same energies (background / efficiency reference) ----
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_4008 = DatasetManager.inclusive_mc.find("703_4009")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4258 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4358 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4416 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# Decay card for the ISR signal, e+e- -> p pbar gamma_ISR, modelled as
# psi(4260) -> p pbar gamma phase space
# (BESIII KKMC convention: psi(4260) as the top mother when no intermediate resonance is produced).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0   p+   anti-p-   gamma     PHSP;
    Enddecay

    End
DECAYCARD

# One 100k-event exclusive MC sample per energy point (identical signal channel, seven datasets).
data_points = [data_3773, data_4008, data_4226, data_4258, data_4358, data_4416, data_4600]
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ppbar_isr_gamma"   # auto-suffixed with the related dataset name
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PpbarISR"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.260]})   # nominal value; ECMS is set per energy point at run time
            .set_alias({"std::vector<double>" => "Vdouble"})

# Build the event selection chain
event_selection = Selection.new
event_selection
  .select_track {                 # charged-track selection: two-prong p pbar topology
      cos_theta   0.93            # |cos(theta)| < 0.93, theta = polar angle of the charged track
      Vz          10.0            # |Vz| < 10 cm
      Vr          1.0             # Vr < 1 cm in the transverse plane
      nChrp       "==1"           # exactly one positively charged track
      nChrn       "==1"           # exactly one negatively charged track
      nNet        "==0"           # net charge zero
  }
  .pid(method: :probability) {    # proton identification with the probability method
      prob_cut   0.001            # PID probability > 0.001 for the tested hypothesis
      identify :proton, against: [:kaon, :pion]   # p and pbar, separated from K and pi
      nprp       "==1"            # exactly one proton
      nprm       "==1"            # exactly one antiproton
  }
  # No explicit photon is required: the ISR photon escapes detection and is inferred from the
  # recoil of the tagged (p, pbar) pair. partial_miss replaces the kinematic fit entirely --
  # this analysis uses no kinematic fit.
  .partial_miss([3]) {            # recID 3 = the undetected photon of the decay card
      # missing mass squared in [-0.1, 0.2] GeV^2/c^4 (sqrt(s) > 4 GeV) -> |m_miss| < 0.447 GeV
      require_recoil_mass 0.0, 0.447
  }

# BOSS-side selection criteria / corrections without a dedicated DSL construct
my_algorithm
  .note(:background_veto,
        "e+ background suppressed by requiring E_EMC/p_rec < 0.5 for the positive track; " \
        "applied after PID and before the missing-momentum / missing-mass cuts")
  .note(:missing_momentum_angle,
        "polar angle of the missing momentum must satisfy theta_miss < 0.125 rad or " \
        "theta_miss > pi - 0.125 rad (undetected ISR photon emitted close to the beam axis)")
  .note(:missing_mass_window,
        "missing mass squared required in [-0.1, 0.2] GeV^2/c^4 for sqrt(s) > 4 GeV and in " \
        "[-0.02, 0.10] GeV^2/c^4 at sqrt(s) = 3.773 GeV; the energy-dependent signed window " \
        "cannot be represented by a single recoil-mass window in the DSL")
  .note(:ppbar_cos_theta_cm,
        "proton and antiproton are additionally required to satisfy |cos theta_{p,pbar}^CM| < 0.75")
  .note(:isr_radiative_correction,
        "NLO Phokhara ISR radiative corrections applied when converting the measured " \
        "e+e- -> p pbar gamma yield into the Born cross section")

# Attach the decay card and render the whole selection into the BOSS C++ algorithm
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on the seven real-data samples, their inclusive MC and the per-energy signal MC
root_files = my_algorithm.execute_on([data_3773, data_4008, data_4226, data_4258,
                                      data_4358, data_4416, data_4600,
                                      incMC_3773, incMC_4008, incMC_4226, incMC_4258,
                                      incMC_4358, incMC_4416, incMC_4600,
                                      *exMCs_signal])