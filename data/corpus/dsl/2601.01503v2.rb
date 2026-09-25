### Dataset description ###
# 4.5 fb^-1 collected at seven c.m. energy points between 4599.53 and 4698.82 MeV
data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

all_data  = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]
all_incMC = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Signal decay card: e+ e- -> Lambda_c+ anti-Lambda_c- with 12 Lambda_c+ decay modes
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000  Lambda_c+  anti-Lambda_c-           PHSP;
    Enddecay

    Decay Lambda_c+
    0.0833  p+   K_S0                          PHSP;
    0.0833  p+   K-   pi+                      PHSP;
    0.0833  p+   K_S0 pi0                      PHSP;
    0.0833  p+   K_S0 pi+  pi-                 PHSP;
    0.0833  p+   K-   pi+  pi0                 PHSP;
    0.0833  Lambda0  pi+                       PHSP;
    0.0833  Lambda0  pi+  pi0                  PHSP;
    0.0833  Lambda0  pi+  pi-  pi+             PHSP;
    0.0833  Sigma0   pi+                       PHSP;
    0.0833  Sigma+   pi0                       PHSP;
    0.0833  Sigma+   pi+  pi-                  PHSP;
    0.0837  p+   pi+  pi-                      PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000  anti-p-  K+   pi-                   PHSP;
    Enddecay

    Decay Sigma0
    1.000  Lambda0  gamma                      PHSP;
    Enddecay

    Decay Sigma+
    0.5157  p+   pi0                           PHSP;
    0.4843  n0   pi+                           PHSP;
    Enddecay

    Decay Lambda0
    1.000  p+  pi-                             PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000  anti-p-  pi+                        PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-                            PHSP;
    Enddecay

    Decay pi0
    1.000  gamma gamma                         PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC generated for each energy point
exMCs = DatasetManager.create_exclusive_mc_for(all_data) do |config, ds|
  config.sample_name     = "LambdacP_hadronic_12modes_#{ds.name}"
  config.related_dataset = ds
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMCs.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) — Lambda_c single tag (12 modes) ###
# NOTE: The paper uses a double-tag technique (both Lambda_c+ and anti-Lambda_c- reconstructed in signal
# channels). DTagTool at BESIII does NOT support Lambda_c double-tag reconstruction, so the DSL below
# expresses the equivalent single-tag reconstruction of anti-Lambda_c- in all twelve hadronic modes.
# The signal-side reconstruction of the recoiling Lambda_c+ is handled in the ROOT analysis stage.
alg_name = "LambdacPHadronicBF"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.682] })
   .note(:analysis_scheme,
         "Double-tag technique used in the paper: both Lambda_c+ (signal) and anti-Lambda_c- (tag) are " \
         "reconstructed in twelve hadronic channels each. BOSS DTagTool does not support Lambda_c double " \
         "tagging, so we generate a single-tag anti-Lambda_c- reconstruction here; the DT combination and " \
         "yields are extracted downstream from the tag NTuple.")
   .note(:intermediate_vetoes,
         "Mass-window vetoes applied to suppress intermediate contributions: (1) M(pbar pi+) in " \
         "[1.17,1.20] GeV/c^2 vetoed in pbar Ks0 pi0, pbar Ks0 pi- pi+, Sigma-bar- pi- pi+, pbar pi- pi+; " \
         "(2) M(pi+ pi-) or M(pi0 pi0) in [0.48,0.52] GeV/c^2 vetoed in Lambda-bar pi- pi+ pi-, " \
         "Sigma-bar- pi0, Sigma-bar- pi- pi+, pbar pi- pi+; " \
         "(3) M(pbar pi0) outside [1.17,1.20] GeV/c^2 required in pbar Ks0 pi0.")
   .note(:best_candidate,
         "For multiple candidates per event, select the one with the minimum |DeltaE|.")

# --- Tag side: anti-Lambda_c- in twelve hadronic modes ---
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,
          :LambdacPtoKPiP,
          :LambdacPtoKsPi0P,
          :LambdacPtoKsPiPiP,
          :LambdacPtoKPiPi0P,
          :LambdacPtoLambdaPi,
          :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi,
          :LambdacPtoSigma0Pi,
          :LambdacPtoSigmaPPi0,
          :LambdacPtoSigmaPPiPi,
          :LambdacPtoPiPiP
  t.charm -1                                       # anti-Lambda_c-
end

# --- Signal side: recoil (all remaining tracks/showers) ---
alg.signal_side do |s|
  s.photons 0..48
end

# --- Kinematic fit: minimal 4C to satisfy the tag_fit requirement ---
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg.apply
root_files = alg.execute_on(all_data + all_incMC + exMCs)
