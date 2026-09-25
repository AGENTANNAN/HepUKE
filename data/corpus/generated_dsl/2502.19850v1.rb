# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
# ψ(2S) (3.686 GeV) real data and inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")
# Off-resonance data at 3.650 GeV, used for the e+e- -> tau+tau- continuum subtraction
cont_data  = DatasetManager.real_data.find("709_3650")

# Decay card: psi(2S) -> tau+ tau-, with tau+ -> e+ nu_e anti-nu_tau and
# tau- -> mu- anti-nu_mu nu_tau.  The tau decays use the EvtGen TAULNUNU model;
# the e+e- -> tau+tau- production (lineshape / ISR) is modelled by ConExc.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0000 ConExc 3;
    Enddecay

    Decay tau+
    1.0000 e+ nu_e anti-nu_tau TAULNUNU;
    Enddecay

    Decay tau-
    1.0000 mu- anti-nu_mu nu_tau TAULNUNU;
    Enddecay

    End
DECAYCARD

# 2M-event exclusive signal MC for psi(2S) -> tau+ tau- (electron + muon final state)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_tautau_emu"
  config.related_dataset = psip_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsiPToTauTau"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {           # exactly two charged tracks: one +, one -
                  cos_theta  0.93        # |cos(theta)| < 0.93
                  Vz         10.0        # |Vz| < 10 cm
                  Vr         1.0         # Vr < 1 cm
                  nChrp      "==1"       # one positively charged track
                  nChrn      "==1"       # one negatively charged track
                  nNet       "==0"       # net charge zero
                }
               .select_photon {          # require zero photons (barrel 25 MeV / endcap 50 MeV)
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam       "==0"       # no photon in the event
                }
               .pid(method: :probability) {   # high-momentum lepton identification
                  # tracks with p > 1.0 GeV/c are treated as leptons; a lepton with
                  # EMC energy > 0.6 GeV is an electron, otherwise a muon
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp  "==1"             # exactly one lepton (l+)
                  nlm  "==1"             # exactly one lepton (l-)
                }
               .kinematic_fit([:lp, :lm]) {   # 4-momentum-constrained fit of the l+ l- system
                  nominal                     # nominal fit
                  constrain_four_momentum
                  chi2_cut 200                # loose chi2 cut; tightened in ROOT
                }

# The e+e- -> tau+tau- continuum is subtracted with the 3.650 GeV data,
# interpolated to the psi(2S) energy; no dedicated DSL construct exists for it.
my_algorithm.note(:continuum_subtraction,
                  "e+e- -> tau+tau- continuum estimated from the 3.650 GeV " \
                  "off-resonance data sample (709_3650), interpolated to the " \
                  "psi(2S) energy (3.686 GeV) via luminosity and the 1/s " \
                  "dependence, and subtracted from the psi(2S) data in ROOT")

# Attach the signal decay card and render the selection into the BOSS algorithm
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Run on psi(2S) data, inclusive MC, the off-resonance continuum data and the signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, cont_data, exMC_signal])