# BESIII: Search for J/psi -> Lambda_c+ e- + c.c.
# Data: 1.31e9 J/psi events at sqrt(s)=3.097 GeV
# Lambda_c+ reconstructed via p K- pi+
# This is a baryon-number-violating decay search.

### Dataset description ###
jpsi_data   = DatasetManager.real_data.find("708_3097")
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for signal J/psi -> Lambda_c+ e- (c.c. implied)
# Lambda_c+ -> p K- pi+
decay_card = <<~DECAYCARD
  Decay J/psi
  1.0000 Lambda_c+ e-                       PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 p+ K- pi+                          PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for signal
exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_Lambdac_e_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

# =============================================================================
# Algorithm: J/psi -> Lambda_c+ e- -> p K- pi+ e-
# Search with 6 competing mass assignment hypotheses (best-chi2 pick at ROOT level)
# =============================================================================
alg = Algorithm.new("JpsiLambdacE")
alg.set_header(["JpsiLambdacEAlg/JpsiLambdacE.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

sel = Selection.new
sel.select_track do
      cos_theta   0.93
      Vz          10.0
      Vr          1.0
      nChrp       "==2"   # exactly 2 positive tracks
      nChrn       "==2"   # exactly 2 negative tracks
      nNet        "==0"   # net charge zero
    end
    # PID: hadrons by per-track probability (highest CL);
    # electron via high-momentum lepton identification
    .pid do
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.1,
                                     treat_as_electron_if_energy_above: 0.6
      identify :proton, against: [:kaon, :pion]
      identify :kaon, against: [:proton, :pion]
      identify :pion, against: [:kaon, :proton]
      nem  "==1"
      nprp "==1"
      nkm  "==1"
      npip "==1"
    end
    # 4C kinematic fit with nominal mass assignment: p K- pi+ e-
    .kinematic_fit([:prp, :km, :pip, :em]) do
      nominal
      constrain_four_momentum
      chi2_cut 200
    end

alg.with_decay_card(decay_card).apply(sel)
alg
  .note(:pid_electron, "Electron: CL_e > 0.001 and CL_e/(CL_e+CL_K+CL_pi) > 0.8. Other tracks assigned by highest CL among pion, kaon, proton hypotheses.")
  .note(:six_mass_hypotheses, "Six mass assignment hypotheses tried: (1) pK-pi+e-, (2) pi+pi-pi+pi-, (3) K+K-K+K-, (4) pi+pi-K+K-, (5) pi+pi-ppbar, (6) K+K-ppbar. Each hypothesis runs its own 4C kinematic fit with the 4 tracks assigned the corresponding masses. Event accepted only if hypothesis (1) gives the best chi2 among all 6. The per-hypothesis mass assignments require separate PID loops in BOSS — the DSL's single-pass PID cannot express this directly, so the competing hypotheses and the best-chi2 veto must be implemented at the ROOT level or in hand-customized BOSS code.")
  .note(:lambda_c_signal_window, "Lambda_c+ mass window: M(pK-pi+) in (2.27, 2.30) GeV/c^2 (~4 sigma resolution). Applied at ROOT level.")
  .note(:result, "No signal events observed. 90% CL upper limit: B(J/psi -> Lambda_c+ e- + c.c.) < 6.9e-8. N_J/psi = 1.31e9. Efficiency = (35.43 +/- 0.02)%. B(Lambda_c+ -> pK-pi+) = (6.35 +/- 0.33)%.")
  .note(:systematic_uncertainties, "Total systematic 7.0%: N_J/psi 0.5%, tracking 1.0%/track, PID (e:0.3%, pi:1.0%, K:0.5%, p:0.6%), kinematic fit 0.2%, B(Lambda_c+) 5.2%. MC modeling negligible.")
  .note(:background, "No background from inclusive J/psi MC. QED background normalized to 0.03 events. Data at off-resonance energies confirm no background in signal region.")
alg.execute_on([jpsi_data, jpsi_incMC, exMC])