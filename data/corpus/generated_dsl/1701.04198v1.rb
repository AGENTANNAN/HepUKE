# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# Real data at the 13 center-of-mass energies from 4.008 to 4.600 GeV
data_points = %w[
  703_4009 703_4090 703_4180 703_4190 703_4200 703_4210 703_4220
  703_4230 703_4260 703_4310 703_4360 703_4420 703_4600
].map { |name| DatasetManager.real_data.find(name) }

# Corresponding inclusive MC for nine of the scan points
incmc_points = %w[
  703_4009 703_4180 703_4190 703_4200 703_4210 703_4220 703_4230
  703_4360 703_4420
].map { |name| DatasetManager.inclusive_mc.find(name) }

# Inclusive Y(4260) background MC at 4.26 GeV (generic psi(4260) decays)
y4260_bkg_incMC = DatasetManager.inclusive_mc.find("703_4260")

# ConExc decay card: e+e- -> p pbar pi0 (mode 50), vhdr -> p+ anti-p- pi0, pi0 -> gamma gamma.
# No `Particle vpho` line for a multi-energy scan: the DSL injects the per-point vpho mass.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1 ConExc 50;
    Enddecay

    Decay vhdr
    1.0000 p+ anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive signal MC at every energy point (ConExc, no KKMC).
# create_exclusive_mc_for shares one card/cross section and derives one sample per energy point.
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ppbarpi0_conexc_exmc"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name  = "ppbarpi0"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})

event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93
    Vz        10.0      # |Vz| < 10 cm
    Vr        1.0       # Vr < 1 cm
    nChrp     "==1"     # exactly one positive charged track
    nChrn     "==1"     # exactly one negative charged track
    nNet      "==0"     # net charge zero
  }
  .select_photon {
    tdc_emc_start     0       # EMC timing window start
    tdc_emc_end       14      # EMC timing window end
    energyThreshold_b 0.025   # > 25 MeV in the barrel
    energyThreshold_e 0.050   # > 50 MeV in the endcap
    nGam              ">=2"   # at least two photons
  }
  .pid(method: :probability) {
    prob_cut 0.001                       # probability cut
    identify :proton, against: [:kaon, :pion]  # separate p+/p- from K and pi
    nprp "==1"                           # exactly one proton
    nprm "==1"                           # exactly one antiproton
  }
  .select_isolated_photon {
    angle_to_prp_track 10.0   # > 10 degrees from the proton track
    angle_to_prm_track 30.0   # > 30 degrees from the antiproton track
    nGam               ">=2"  # still at least two photons
  }
  # Nominal 4C kinematic fit to p pbar gamma gamma, constraining the total four-momentum
  # to the e+e- CMS. Loose BOSS chi2 cut; the tighter chi2 < 30 and the pi0 mass window
  # |M(gg) - m_pi0| < 15 MeV/c^2 are applied downstream in ROOT.
  .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

# Generate the BOSS algorithm for the signal process and run over all datasets
algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = algorithm.execute_on(data_points + incmc_points + [y4260_bkg_incMC] + exMC_signal)