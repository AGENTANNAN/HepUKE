# Paper: 2304.14655v1
# Title: Test of CP symmetry using entangled Sigma+ anti-Sigma- pairs from J/psi decays
# Energy: 3.097 GeV (J/psi, single energy point)
# Final state: J/psi -> Sigma+ anti-Sigma-; Sigma+ -> n pi+ or p pi0;
#             anti-Sigma- -> nbar pi- or pbar pi0
# Measurement of CP-odd weak decay parameters alpha+ and alpha_bar-
# Neutron/anti-neutron reconstructed via EMC showers
# Similar to 2209.14564v2 but with Sigma instead of psi(3686)

### Dataset preparation ###
data_708_3097 = DatasetManager.load_real_data.find("708_3097")
incMC_708_3097 = DatasetManager.load_inclusive_mc.find("708_3097")

all_data = [data_708_3097]
all_incMC = [incMC_708_3097]

# Decay card: J/psi -> Sigma+ Sigma- bar
# Sigma+ -> n pi+; Sigma- -> nbar pi-
# Note: neutron/anti-neutron reconstructed via EMC showers
decay_card = <<~DECAYCARD
    Decay J/psi
    1.000  Sigma+  anti-Sigma-                 HELAMP 1.0 0.0;
    Enddecay

    Decay Sigma+
    1.000  n0  pi+                              PHSP;
    Enddecay

    Decay anti-Sigma-
    1.000  anti-n0  pi-                         PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Jpsi_Sigma_Sigmabar_CP"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("Jpsi_SigmaSigmabar_CP")
alg.set_header(["JpsiSigmaCPAlg/JpsiSigmaSigmabarCP.h"])
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new

# Charged tracks: exactly 2 (pi+ and pi-), zero net charge
# Neutron/anti-neutron reconstructed via EMC (not charged tracks)
event_selection.select_track {
                 cos_theta 0.93
                 Vz   30.0         # |Vz| < 30 cm
                 Vr   10.0         # |Vxy| < 10 cm
                 nTot "==2"        # 2 charged pions, zero net charge
               }
               # PID: both tracks identified as pions
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :pion, against: [:kaon, :proton]
                 npip "==1"
                 npim "==1"
               }
               .assign({:chrgp => :pip, :chrgn => :pim})
               # 4C kinematic fit
               .kinematic_fit([:pip, :pim]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 200
               }

alg.with_decay_card(decay_card).apply(event_selection)

# Neutron/anti-neutron reconstruction via EMC showers
# - EMC shower energy > 600 MeV (barrel/endcap)
# - Lateral moment > 20 (photon suppression)
# - Angle to charged track > 10 deg
# - EMC time [0, 700] ns
alg.note(:neutron_emc_selection,
  "Neutron/anti-neutron reconstructed via EMC showers: E(EMC) > 600 MeV; lateral moment > 20; angle to charged track > 10 deg; EMC time in [0,700] ns. Most energetic candidate selected. Applied in ROOT.")

# Kinematic fit with free neutron parameters
# Four-momentum conservation + Sigma mass constraints
alg.note(:kinematic_fit_free_neutron,
  "Kinematic fit: 4-momentum conservation + Sigma+ and Sigma- mass constraints. Neutron momenta free parameters; anti-neutron angles used, energy floated. Applied in ROOT.")

# Entangled Sigma+ anti-Sigma- pairs from J/psi for CP test
# Measure decay parameters alpha+ (Sigma+ -> n pi+) and alpha_bar- (Sigma- -> nbar pi-)
# CP asymmetry: A_CP = (alpha+ + alpha_bar-) / (alpha+ - alpha_bar-)
alg.note(:cp_measurement,
  "CP test using entangled Sigma+ anti-Sigma- pairs. alpha+ = -0.0565+/-0.0047_stat+/-0.0022_syst; alpha_bar- = 0.0481+/-0.0031_stat+/-0.0019_syst. A_CP = -0.080+/-0.052_stat+/-0.028_syst. First measurement of alpha- decay parameter. Applied in ROOT using angular analysis.")

# Also uses Sigma+ -> p pi0 and Sigma- -> pbar pi0 channels
# for decay parameter ratios relative to reference modes
alg.note(:additional_channels,
  "Also uses Sigma+ -> p pi0 and Sigma- -> pbar pi0 decays for ratio measurements: alpha+/alpha0 and alpha_bar-/alpha_bar0. Pi0 reconstructed from gamma gamma with mass-constrained kinematic fit. Applied in ROOT.")

# 1.0087 x 10^10 J/psi events
alg.note(:total_jpsi_events,
  "Total of (1.0087 +/- 0.0044) x 10^10 J/psi events used.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)