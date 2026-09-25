# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# High-luminosity XYZ scan points (sqrt(s) = 3.773 - 4.600 GeV) for the
# e+e- -> pi+pi-J/psi Born cross-section measurement.
data_3773 = DatasetManager.real_data.find("712_3773")   # psi(3770)
data_3810 = DatasetManager.real_data.find("703_3810")
data_3900 = DatasetManager.real_data.find("703_3900")
data_4009 = DatasetManager.real_data.find("703_4009")
data_4090 = DatasetManager.real_data.find("703_4090")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4360 = DatasetManager.real_data.find("703_4360")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4470 = DatasetManager.real_data.find("703_4470")
data_4600 = DatasetManager.real_data.find("703_4600")

data_points = [data_3773, data_3810, data_3900, data_4009, data_4090, data_4190,
               data_4210, data_4220, data_4230, data_4260, data_4360, data_4420,
               data_4470, data_4600]

# Corresponding inclusive MC samples at the same energy points
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_3810 = DatasetManager.inclusive_mc.find("703_3810")
incMC_3900 = DatasetManager.inclusive_mc.find("703_3900")
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")
incMC_4090 = DatasetManager.inclusive_mc.find("703_4090")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4470 = DatasetManager.inclusive_mc.find("703_4470")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

incMC_samples = [incMC_3773, incMC_3810, incMC_3900, incMC_4009, incMC_4090,
                 incMC_4190, incMC_4210, incMC_4220, incMC_4230, incMC_4260,
                 incMC_4360, incMC_4420, incMC_4470, incMC_4600]

# ConExc decay cards (mode 90) for the continuum process e+e- -> pi+pi-J/psi.
# ConExc models ISR up to second order; the DSL auto-detects the "ConExc" token,
# switches to the no-KKMC simulation template and injects "Particle vpho <ECMS> 0.0"
# per energy point (no psi(4260) KKMC mother, no explicit Particle vpho line).
# Two cards cover the two J/psi -> l+l- signal channels (mu+mu- and e+e-).
decay_card_mumu = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 90;
    Enddecay

    Decay J/psi
    1.0 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_ee = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 90;
    Enddecay

    Decay J/psi
    1.0 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Generate 60k-event exclusive signal MC at every scan point. The card / cross
# section / event count are shared; only the related dataset changes, so the
# batch creator runs one MC per energy point for each lepton channel.
exMCs_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_pipiJpsi_mumu"
  config.events      = 60_000
  config.decay_card  = decay_card_mumu
  config.cross_section = :default
end

exMCs_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_pipiJpsi_ee"
  config.events      = 60_000
  config.decay_card  = decay_card_ee
  config.cross_section = :default
end

# Persist the MC configurations (one file per energy point per channel)
exMCs_mumu.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }
exMCs_ee.each   { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
alg_name = "pipiJpsi"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.26]})   # per-point ECMS injected by ConExc; 4.26 as nominal value
   .set_alias({"std::vector<double>" => "Vdouble"})

# Single selection chain shared by both J/psi -> l+l- channels (leptons handled
# generically through the lp / lm lists).
event_selection = Selection.new
  .select_track {                 # Charged track quality + multiplicity
      cos_theta 0.93              # |cos(theta)| < 0.93
      Vz        10.0              # |Vz| < 10 cm (beam direction)
      Vr        1.0               # Vr < 1 cm (transverse plane)
      nChrp     "==2"             # exactly 2 positive tracks
      nChrn     "==2"             # exactly 2 negative tracks (==> exactly four charged tracks)
      nNet      "==0"             # net charge zero
  }
  .select_photon {                # Photon selection (photons are stored but not put in the fit)
      tdc_emc_start     0         # EMC TDC window 0-14
      tdc_emc_end       14
      energyThreshold_b 0.025     # 25 MeV barrel threshold
      energyThreshold_e 0.050     # 50 MeV endcap threshold
      angle_to_track    10.0      # >10 degrees from nearest charged track
  }
  .pid(method: :probability) {    # Probability PID, prob cut 0.001
      prob_cut 0.001
      # Treat tracks with p > 1.06 GeV/c as leptons; a lepton is called electron
      # if its EMC energy is above 0.6 GeV, otherwise muon.
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.06,
                                     treat_as_electron_if_energy_above: 0.6
      identify :pion, against: [:kaon]   # pi/K separation (pi+ and pi-)
      npip      "==1"
      npim      "==1"
      nlp       "==1"
      nlm       "==1"
  }
  # 4C kinematic fit to pi+pi-l+l- (constrain the final-state four-momentum to CMS)
  .kinematic_fit([:pip, :pim, :lp, :lm]) {
      nominal                     # nominal fit: fitted four-momenta saved
      constrain_four_momentum
      chi2_cut 60
  }

# Attach the decay card and render the event-selection algorithm
alg.with_decay_card(decay_card_mumu).apply(event_selection)

# Execute on real data, inclusive MC and the per-point exclusive signal MCs
root_files = alg.execute_on(data_points + incMC_samples + exMCs_mumu + exMCs_ee)