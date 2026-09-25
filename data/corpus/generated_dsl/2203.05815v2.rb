### Dataset description ###
# Energy scan 4.2263-4.6984 GeV:  e+e- -> pi+ pi- psi2(3823)
#   psi2(3823) -> gamma chi_c1, chi_c1 -> gamma J/psi, J/psi -> l+ l- (e+e-  or mu+mu-)
# 20 real-data points and their matching inclusive-MC samples
scan_points = %w[
  703_4230 703_4237 703_4246 703_4260 703_4270 703_4280 703_4310 703_4360
  703_4390 703_4420 703_4470 703_4530 703_4575 703_4600
  706_4612 706_4620 706_4640 706_4660 706_4680 706_4700
]
data_points  = scan_points.map { |name| DatasetManager.real_data.find(name) }
incMC_points = scan_points.map { |name| DatasetManager.inclusive_mc.find(name) }

# Decay cards (EvtGen) for the full decay chain, one per J/psi lepton channel
decay_card_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- psi(3823) PHSP;
  Enddecay

  Decay psi(3823)
  1.0000 gamma chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.0000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

decay_card_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- psi(3823) PHSP;
  Enddecay

  Decay psi(3823)
  1.0000 gamma chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.0000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# 50k-event exclusive MC at every scan point (one MC per energy point), both lepton channels
exMCs_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_psi3823_scan_ee"    # auto-suffixed per energy point
  config.events        = 50_000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

exMCs_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_psi3823_scan_mumu"  # auto-suffixed per energy point
  config.events        = 50_000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "Psi3823ToGammaChiC1"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.23]})   # representative scan energy (per-point value set at execute_on)
            .set_alias({"std::vector<double>" => "Vdouble"})
            # track helix-parameter corrections applied to MC before the kinematic fit
            .note(:helix_correction, "track helix-parameter corrections applied to simulated (MC) tracks before the 4C kinematic fit")
            # gamma-conversion veto (angular criterion, not expressible as a mass window)
            .note(:background_veto, "gamma-conversion veto: reject events with cos(opening angle between pi+ and pi-) < 0.98")

# Common selection for both lepton channels (lepton universality -> combined l+/l- lists)
event_selection = Selection.new
  .select_track {                    # exactly two positive and two negative charged tracks
    cos_theta 0.93                   # |cos(theta)| < 0.93
    Vz        10.0                   # |Vz| < 10 cm
    Vr        1.0                    # Vr < 1 cm
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {                   # two signal photons (psi2 -> gamma chi_c1 -> gamma J/psi)
    tdc_emc_start 0
    tdc_emc_end   14
    angle_to_track 10.0
    energyThreshold_b 0.025          # 25 MeV in the barrel
    energyThreshold_e 0.050          # 50 MeV in the endcap
    nGam ">=2"
  }
  .pid(method: :probability) {       # PID by the probability method
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6  # e+- and mu+- combined as l+/l-
    identify :pion, against: [:kaon]  # pi+ / pi- separated from kaons
    npip "==1"
    npim "==1"
    nlp  "==1"
    nlm  "==1"
  }
  .kinematic_fit([:pip, :pim, :lp, :lm, :gamma, :gamma]) {  # fit to pi+ pi- l+ l- gamma gamma
    nominal                        # nominal fit -> corrected four-momenta are stored
    constrain_four_momentum        # 4C constraint on the CMS four-momentum
    chi2_cut 200
    invariant_mass_of(:lp, :lm).within(3.06, 3.135)                  # J/psi mass window
    invariant_mass_of(:pip, :pim, :lp, :lm).out_of(3.679, 3.693)     # psi(2S) veto: |M(pi+pi-l+l-) - m[psi(2S)]| > 7 MeV
    invariant_mass_of(:gamma, :gamma, :pip, :pim).larger_than(0.65)  # eta J/psi veto
  }

# One algorithm shared by both lepton channels (identical final-state topology and selection)
my_algorithm.with_decay_card(decay_card_ee).apply(event_selection)

# Execute on all real-data points, their inclusive MC, and the two exclusive-MC batches
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs_ee + exMCs_mumu)