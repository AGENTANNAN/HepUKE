# Amplitude analysis of the SCS decay Lambda_c+ -> p K+ K-
# 4.4 fb^-1 e+e- data at CM energies 4600-4698 MeV
# Single-tag Lambda_c+ reconstructed directly via p K+ K- final state,
# identified by M_BC and Delta_E variables.

### Datasets (approximate matches for 4600, 4628, 4641, 4661, 4682, 4698 MeV) ###
data_4600 = DatasetManager.real_data.find("703_4600")
data_4620 = DatasetManager.real_data.find("706_4620")   # 4628 MeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4641 MeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4661 MeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4682 MeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4698 MeV

data_points = [data_4600, data_4620, data_4640, data_4660, data_4680, data_4700]

incMC_points = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
]

# Signal MC uses PHSP for amplitude analysis input; final signal MC uses amplitude model
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 Lambda_c+ anti-Lambda_c-               PHSP;
    Enddecay

    Decay Lambda_c+
    1.000 p+  K+  K-                             PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000 anti-p-  K-  K+                        PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name    = "Lambdac_pKK_signal_PHSP"
  config.events          = 1000000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: "exMC_LambdacpKK_#{m.name}") }

alg = Algorithm.new("LambdacpKK")
alg.set_header(["LambdacpKKAlg/LambdacpKK.h"])
   .set_constant({"ECMS" => [:double, 4.682]})   # nominal; per-run energy handled at runtime

sel = Selection.new
sel.select_track {
       cos_theta 0.93       # |cos(theta)| < 0.93
       Vr        1.0        # |Vr| < 1 cm for K
       Vz        10.0       # |Vz| < 10 cm
       nChrp     ">=2"      # p (or anti-p) + K+ (or K-)
       nChrn     ">=1"
     }
    .pid(method: :probability) {
       prob_cut 0.0
       identify :proton, against: [:kaon, :pion]
       identify :kaon,   against: [:pion, :proton]
       identify :pion,   against: [:kaon, :proton]
       nprp ">=1"
       nkp  ">=1"
       nkm  ">=1"
     }
    # Additional tight vertex cut on proton: Vr < 0.2 cm handled in ROOT
    # Kinematic fit: constrain to 4-momentum of e+ e- CM
    .kinematic_fit([:prp, :kp, :km]) {
       nominal
       constrain_four_momentum
       chi2_cut 200
     }

alg
  .note(:proton_tight_Vr,
        "Additional |Vr| < 0.2 cm cut applied to proton candidate to suppress secondary/beam-related protons.")
  .note(:MBC_preselection,
        "M_BC > 2.25 GeV/c^2 and |Delta_E| < 0.1 GeV loose preselection; Delta_E = E - E_beam, M_BC = sqrt(E_beam^2 - |p|^2).")
  .note(:best_candidate_selection,
        "If multiple pK+K- combinations pass, retain the one with the smallest |Delta_E|.")
  .note(:deltaE_windows,
        "For amplitude analysis: |Delta_E| < 0.005 GeV. For branching-fraction fit: |Delta_E| < 0.02 GeV.")
  .note(:sPlot_signal_extraction,
        "MBC fit uses MC-simulated signal shape convolved with Gaussian; background ARGUS function with endpoint = E_beam; sWeights via sPlot for amplitude analysis.")
  .note(:amplitude_components,
        "Resonant amplitude components: phi(1020) [Breit-Wigner, mass/width free], Lambda(1670) [BW with PDG values], f0(980) [Flatte pi pi/KK], Lambda(1405) [Flatte pK/Sigma pi].")
  .note(:multi_energy_simultaneous_fit,
        "Simultaneous MBC fit at six CM energies (4600, 4628, 4641, 4661, 4682, 4698 MeV) with common branching fraction; Gaussian resolution parameters shared.")
  .with_decay_card(decay_card_signal)
  .apply(sel)
alg.execute_on(data_points + incMC_points + exMC_signal)
