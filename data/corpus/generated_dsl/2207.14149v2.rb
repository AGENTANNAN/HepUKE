# =============================================================================
# Λc+ → Λ e+ νe  (double-tag: 14 hadronic Λc− tag modes + signal side)
# Seven c.m. energy points from 4.600 to 4.699 GeV (total 4.5 fb⁻¹)
# Tag-based analysis → TagAnalysis (no Selection object)
# =============================================================================

### ------------------------------ Datasets ----------------------------------
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.5995 GeV
data_4610 = DatasetManager.real_data.find("706_4610")   # 4.6119 GeV
data_4620 = DatasetManager.real_data.find("706_4620")   # 4.6280 GeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4.6409 GeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4.6612 GeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4.6819 GeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4.6988 GeV
data_points = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]

incMC_points = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700")
]

### --------------------- Decay card (signal mode) ---------------------------
# Λc+ → Λ e+ νe,  Λ → p π−   (top mother = psi(4260), BESIII KKMC convention)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c-   PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Lambda0 e+ nu_e            PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                     HypWK;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K+ pi-             PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events, one sample per energy point (shared card/events)
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_LambdacToLambdaEnu"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### ------------------------- Event selection (BOSS) -------------------------
alg_name = "LambdacToLambdaEnu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   # nominal energy of the 4.600–4.699 GeV scan; the tag fit takes the per-run
   # measured beam energy (dtag_reconstruction beam_energy :db default)
   .set_constant({ "ECMS" => [:double, 4.640] })
   .with_decay_card(decay_card_signal)
   .note(:lambda_vertex_fit,
         "signal-side Λ → p π− is required to form a common secondary vertex; the tag " \
         "layer exposes no vertex-fit primitive, so the Λ is selected by the " \
         "|M(p π−) − M(Λ)| < 5 MeV mass window applied in the kinematic fit instead")
   .note(:electron_pid_criteria,
         "e+ identified with CL_e > 0.001, CL_e > CL_π and CL_e > CL_K; the signal-side " \
         "electron uses the fixed v1 lepton PID recipe of the tag layer (not DSL-tunable)")

# --- Tag side: Λc− reconstructed from the 14 hadronic tag modes (ST) ---------
alg.tag_side(:Lambdac) do |t|
  # the analysed mode set is the hadronic Λc− tag group; the modes explicitly
  # quoted in the analysis are pK−π+, Λπ−, Λπ−π+π− and pK−π+π−π+
  t.mode_group :hadronic
  t.charm(-1)                      # pin the tagged Λc− side
end

# --- Signal side: what the tag did not use ----------------------------------
# one p and one π− (Λ → p π−) plus one e+; the neutrino escapes undetected
alg.signal_side do |s|
  s.charged(prp: 1, pim: 1, ep: 1)
  s.missing :nu_e                  # massless missing particle (semileptonic νe)
end

# --- Kinematic fit: 4C constraint over tag + signal + missing ---------------
alg.fit do |f|
  f.constrain_four_momentum                    # 4-momentum conservation to the measured CMS 4-vector
  f.invariant_mass_of(:prp, :pim).between(1.1107, 1.1207)  # |M(p π−) − M(Λ)| < 5 MeV
  f.chi2_cut 200                               # loose χ² cut; optimal cut applied in ROOT
  # U_miss = E_miss − |p_miss| (and q² for the single signal lepton) are stored
  # automatically with the missing-particle declaration; the −0.06 < U_miss < 0.06 GeV
  # window, the 4D ML fit and the five form-factor parameters are ROOT-level steps.
end

alg.apply                                  # tag spec: apply takes no Selection argument
root_files = alg.execute_on(data_points + incMC_points + exMC_signal)