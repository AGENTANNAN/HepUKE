# Paper: arXiv:2401.14711v1
# e+e- -> pi+ pi- pi0 at 19 CM energies (2.00-3.08 GeV)
# Ordinary analysis with ConExc generator and 4C kinematic fit

### Dataset preparation ###

# Energy scan from 2.00 to 3.08 GeV (19 points) - R-scan data
# Representative query for R-scan data points
rscan_data = DatasetManager.real_data.where({ "Ecms(MeV)" => (2000..3080) })

# ConExc decay card for e+e- -> pi+ pi- pi0
# Mode 7 = pi+ pi- pi0 (1.0625-2.9875 GeV, no phi, no J/psi)
conexc_decay_card = <<~DECAYCARD
    Decay vpho
    1 ConExc 7;
    Enddecay
    Decay vhdr
    1 pi+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC for all scan points (omit Particle vpho for per-point injection)
sig_mc_scan = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name = "exmc_pipipi0_scan"
  config.events = 50000
  config.decay_card = conexc_decay_card
  config.cross_section = :default
end

### Event selection ###

alg = Algorithm.new("PiPiPi0Scan")
alg.set_header(["PiPiPi0ScanAlg/PiPiPi0Scan.h"])
   .set_constant({"ECMS" => [:double, 2.50]})   # placeholder; multi-energy analysis uses per-run beam energy

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=1"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
    npim ">=1"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:pip, :pim, :pi0]) do
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

# Notes
alg.note(:conexc_mode, "ConExc mode 7 used for pi+ pi- pi0 continuum production (includes Born + ISR up to 2nd order)")
alg.note(:energy_scan, "Analysis at 19 CM energies from 2.00 to 3.08 GeV; ECMS set per-run via MeasuredEcmsSvc")
alg.note(:track_cuts, "Additional track cuts: E/p < 0.8 for electron rejection; opening angle theta(pi+,pi-) < 160 degrees (applied in post-fit analysis)")
alg.note(:pi0_signal_region, "pi0 signal region |M(gamma gamma) - M_pi0| < 15 MeV/c^2 (applied in ROOT, not in DSL)")
alg.note(:helicity_cut, "Helicity angle |cos(theta_h)| < 0.8 cut applied in PWA stage (ROOT analysis, beyond DSL scope)")
alg.note(:pwa, "Partial wave analysis (PWA) performed on accepted events; yield extracted from M(gamma gamma) fit (ROOT analysis, beyond DSL scope)")
alg.note(:born_cross_section, "Born cross section measurement using ISR correction factor f_ISR from ConExc generator log")

alg.with_decay_card(conexc_decay_card).apply(event_selection)
alg.execute_on([rscan_data, sig_mc_scan].flatten)