### Dataset description ###
# Continuum data at sqrt(s) = 3.65 GeV (62 pb^-1) and the matching inclusive MC
data_3650  = DatasetManager.real_data.find("709_3650")
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")

# --- Continuum MC decay card (ConExc, "ConExc/psi(4260)" convention) -------
# e+e- -> qqbar (u, d, s) continuum generated with the ConExc model.  The literal
# token `ConExc` makes the DSL switch to the no-KKMC simulation template and
# inject `Particle vpho <ECMS> 0.0`, so `Particle vpho` is deliberately omitted.
decay_card_continuum = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 1;
    Enddecay
    End
DECAYCARD

# 500k-event exclusive continuum MC for e+e- -> qqbar (u, d, s)
exMC_continuum = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3650_qqbar_continuum"
  config.related_dataset = data_3650
  config.events          = 500000
  config.decay_card      = decay_card_continuum
  config.cross_section   = :default
end

# --- Signal decay card for the algorithm header ---------------------------
# Reconstructed signal final state: the charged-dipion system produced in the
# continuum; the unmeasured remainder X is not modelled.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi-  PHSP;
    Enddecay
    End
DECAYCARD

### Event selection (BOSS) ###
alg_name = "PiPiX"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.65]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
      cos_theta 0.93        # |cos(theta)| < 0.93 for charged tracks
      Vz        10.0        # |Vz| < 10 cm
      Vr        1.0         # Vr < 1 cm
      nTot      ">=3"       # at least three charged tracks in total
    }
   .select_photon {
      tdc_emc_start     0   # cluster timing within 0-14 TDC counts (delay <= 700 ns)
      tdc_emc_end       14
      energyThreshold_b 0.025  # 25 MeV in the barrel
      energyThreshold_e 0.050  # 50 MeV in the endcup
    }
   .pid(method: :probability) {
      prob_cut 0.001        # PID probability cut
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6  # p>1.0 -> lepton; EMC raw E>0.6 -> electron
      identify :pion, against: [:kaon, :proton]  # pi-K-p separation
      npip ">=1"            # at least one pi+
      npim ">=1"            # at least one pi-
    }

# Selection criteria / procedures that the current DSL cannot express
my_algorithm
  .note(:inclusive_no_kinematic_fit, "the reaction is semi-inclusive, e+e- -> pi+pi- X; the undetected system X is not reconstructed, so no kinematic fit is applied")
  .note(:electron_veto, "electron veto on charged tracks: require L(e) > 0.001 and L(e)/(L(e)+L(pi)+L(K)) > 0.8 using the PID likelihoods")
  .note(:visible_energy, "require the total visible energy (sum over all selected charged tracks and photons) > 1.5 GeV")
  .note(:pion_fractional_energy, "require both pion scaled energies zi = 2 E_pi / sqrt(s) to lie within [0.2, 0.9]")
  .note(:opening_angle, "require the opening angle between the two pion candidates > 120 degrees")
  .note(:dipion_pairing, "if more than two pions survive, form all pi+pi- combinations; the two pions of each pair are randomly labelled h1/h2 and a single pion may enter several pairs (evaluated in the ROOT analysis)")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([data_3650, incMC_3650, exMC_continuum])