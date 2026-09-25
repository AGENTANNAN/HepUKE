# Paper: arXiv:2401.09468v2
# e+e- -> Sigma+ Sigma- at 41 CM energies (3.510-4.951 GeV)
# Sigma+ -> p pi0, Sigma- -> anti-p pi0, pi0 -> gamma gamma
# Ordinary analysis with ConExc generator and 6C kinematic fit

### Dataset preparation ###

# Load R-scan / scan data points at 41 energies
# Representative energy points (full list of 41 would be loaded via group_by_boss query)
# Energy scan from 3.510 to 4.951 GeV
data_points = DatasetManager.real_data.where({ "Ecms(MeV)" => (3510..4951) })

incMC_by_energy = {}  # Inclusive MC per scan point (simplified)

# ConExc decay card for e+e- -> Sigma+ anti-Sigma+
# Mode 44 = Sigma pair production (charged Sigma pairs, 2.308-5.000 GeV)
# vhdr daughter names: Sigma+ and anti-Sigma+ (EvtGen names)
conexc_decay_card = <<~DECAYCARD
    Decay vpho
    1 ConExc 44;
    Enddecay
    Decay vhdr
    1 Sigma+ anti-Sigma+ PHSP;
    Enddecay
    Decay Sigma+
    1 p+ pi0 PHSP;
    Enddecay
    Decay anti-Sigma+
    1 anti-p- pi0 PHSP;
    Enddecay
    Decay pi0
    1 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC: one card for all scan points (omit Particle vpho)
# Using create_exclusive_mc_for for multi-energy scan
sig_mc_scan = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_SigmaSigma_scan"
  config.events = 50000
  config.decay_card = conexc_decay_card
  config.cross_section = :default
end

### Event selection ###

alg = Algorithm.new("SigmaSigmaScan")
alg.set_header(["SigmaSigmaScanAlg/SigmaSigmaScan.h"])
   .set_constant({"ECMS" => [:double, 4.26]})   # placeholder; multi-energy analysis uses per-run beam energy

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"    # proton + pi+ from Sigma+ decay... wait, actually just p+ and pbar from Sigma decays
    nChrn ">=2"    # anti-p-
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=4"     # two pi0 -> 4 photons
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  end
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

# Notes
alg.note(:conexc_mode, "ConExc mode 44 used for Sigma+ anti-Sigma+ pair production (charged Sigma pairs)")
alg.note(:energy_scan, "Analysis at 41 CM energies from 3.510 to 4.951 GeV; ECMS set per-run via MeasuredEcmsSvc")
alg.note(:proton_pid, "Proton ID: momentum > 0.5 GeV/c (applied as additional cut beyond DSL PID block)")
alg.note(:signal_sideband, "Signal/sideband method used in ROOT for yield extraction with iterative weighting for ISR correction")
alg.note(:iterative_isr, "Iterative weighting procedure applied to correct for ISR effects (ROOT analysis)")
alg.note(:born_cross_section, "Born cross section measurement using ISR correction factor f_ISR from ConExc generator log")
alg.note(:custom_xs, "ConExc cross section from internal Born table; f_vacuum from generator log used in sigma calculation")

alg.with_decay_card(conexc_decay_card).apply(event_selection)
alg.execute_on([data_points, sig_mc_scan].flatten)