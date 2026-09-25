# Analysis: Measurement of Born cross sections for e+ e- -> Sigma- Sigmabar+
# at sqrt(s) = 3.51-4.95 GeV, and observation of psi(3770) -> Sigma- Sigmabar+.
# Reconstruction: Sigma- -> n pi- (neutron NOT reconstructed);
#                 Sigmabar+ -> nbar pi+ (only the antineutron detected in the EMC).
# Partial-reconstruction technique with the neutron treated as missing.

### Dataset description ###
# Representative selection of c.m. energy points covered by BESIII from
# sqrt(s) = 3.51 to 4.95 GeV (the paper uses 52 points; only the primary
# psi(3770) and the XYZ/scan energies commonly available in DatasetManager
# are listed here).
scan_datasets = [
  DatasetManager.real_data.find("712_3773"),
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4914"),
  DatasetManager.real_data.find("707_4946"),
]

scan_incMCs = [
  DatasetManager.inclusive_mc.find("712_3773"),
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4914"),
  DatasetManager.inclusive_mc.find("707_4946"),
]

# ConExc decay card for e+e- -> Sigma- Sigmabar+ (mode 44).
# 'Particle vpho' is intentionally omitted so that the DSL auto-injects the
# correct sqrt(s) at every scan energy.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1 ConExc 44;
    Enddecay

    Decay vhdr
    1 Sigma- anti-Sigma+                      PHSP;
    Enddecay

    Decay Sigma-
    1.0000 n0 pi-                             PHSP;
    Enddecay

    Decay anti-Sigma+
    1.0000 anti-n0 pi+                        PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC at every scan point (100 000 events per point in paper).
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name   = "SigmaMSigmaP_ConExc"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default   # inert for ConExc but required
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
alg_name = "SigmaMSigmaP"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 3.773] })   # nominal point
            .set_alias({ "std::vector<double>" => "Vdouble" })

# Selection chain: one pi+ and one pi- track; one antineutron EMC shower.
event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93          # |cos(theta)| < 0.93
                  Vz        30.0          # |Vz|  < 30 cm  (loose IP cut)
                  Vr        10.0          # |Vxy| < 10 cm  (loose IP cut)
                  nChrp    "==1"
                  nChrn    "==1"
                  nNet     "==0"
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14    # EMC time within [0, 700] ns
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  angle_to_track    15.0
                  nGam ">=1"
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :pion, against: [:kaon, :proton]
                  npip "==1"
                  npim "==1"
                }
               # Partial reconstruction: reconstruct Sigmabar+ from (nbar, pi+)
               # while the neutron on the Sigma- side is inferred from recoil.
               .partial_miss([1]) {         # miss the Sigma- side
                  require_recoil_mass 1.10, 1.35
                }

my_Algorithm
  .note(:antineutron_identification,
        "Antineutron identified as the EMC shower with the highest deposited " \
        "energy; deposited energy > 0.62 GeV at sqrt(s) = 3.773 GeV (varies " \
        "per energy point). Shower second moment S = sum E_i r_i^2 / sum E_i " \
        "> 18 cm^2; total number of hits associated with the shower > 16; " \
        "angle between the shower and any charged track > 15 degrees.")
  .note(:cosmic_veto,
        "Difference of pi+ and pi- times of flight in the MDC required to be " \
        "greater than -5 ns to remove cosmic-ray background.")
  .note(:one_c_kinematic_fit,
        "1-C kinematic fit with the (pi+ nbar) invariant mass constrained to " \
        "the Sigmabar+ nominal mass and with the antineutron flight direction " \
        "constrained by the EMC shower position; chi2_1C < 60.  The recoil " \
        "mass M(pi+ nbar) is required to lie in [1.10, 1.35] GeV/c^2.")
  .note(:iterative_efficiency,
        "Detection efficiency and ISR correction factor determined via an " \
        "iterative procedure that inputs the measured Born cross section.")

my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_Algorithm.execute_on(scan_datasets + scan_incMCs + exMC_signal)
