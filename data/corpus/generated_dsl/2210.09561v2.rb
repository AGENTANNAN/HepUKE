# ============================================================================
# Inclusive Λc− → n̄ X measurement
# Single-tag method: tag Λc+ → p K− π+
# Seven c.m. energy points: 4.600 – 4.700 GeV
# ============================================================================

### Dataset preparation ###
# Real data at the seven energy points (BOSS version _ CMS energy)
data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
real_data = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]

# Matching inclusive MC samples for the seven energy points
incMC = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700")
]

# Decay card for the signal process (EvtGen format, EvtGen particle names).
# e+e- → Λc+ Λc− ; tag Λc+ → p K− π+ ; signal Λc− → n̄ X (inclusive, PHSP proxy).
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0 p+ K- pi+ PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0 anti-n- pi- PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for the Λc+Λc− background: the opposite Λc− decays to p̄ X and is
# misidentified as an antineutron (PHSP proxy for the inclusive p̄ X final state).
decay_card_bkg = <<~DECAYCARD
  Decay psi(4260)
  1.0 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0 p+ K- pi+ PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0 anti-p- pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Signal exclusive MC — same signal MC generated at each of the seven energy points
exMCs_signal = DatasetManager.create_exclusive_mc_for(real_data) do |config|
  config.sample_name   = "exmc_lambdac_nbarX"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# Λc+Λc− background exclusive MC at each energy point
exMCs_bkg = DatasetManager.create_exclusive_mc_for(real_data) do |config|
  config.sample_name   = "exmc_lambdacp_lambdacm_pbarX"
  config.events        = 100000
  config.decay_card    = decay_card_bkg
  config.cross_section = :default
end

### Tag-based event selection (BOSS) ###
# Single tag: reconstruct Λc+ from the pre-stored DTag collection, pick out the
# antineutron recoiling against it; the signal side provides the n̄ and the vetoes.
alg = TagAnalysis.new("LcTagNbarInclusive")
alg.set_header(["LcTagNbarInclusiveAlg/LcTagNbarInclusive.h"])
   .set_constant({ "ECMS" => [:double, 4.600] })
   .with_decay_card(decay_card_signal)

# Tag side: Λc+ → p K− π+
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP                        # Λc+ → p K− π+
  t.charm 1                                      # pin the Λc+ (charm +1) tag
  t.window :deltaE, min: -0.034, max: 0.020      # ΔE window (−34, +20) MeV
  t.window :mBC,    min: 2.275,  max: 2.310      # MBC ∈ (2.275, 2.31) GeV/c²
end

# Signal side: no photons, no protons, antineutron inferred in the fit
alg.signal_side do |s|
  s.photons 0            # require no photons on the signal side
  s.missing :n_bar       # antineutron — momentum floated in the 4C fit
  s.min_photon_angle 20.0 # shower–track opening angle > 20°
end

# 4-constraint kinematic fit
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# Inexpressible BOSS-side procedures preserved as notes
alg
  .note(:proton_veto,
        "signal side: no charged track may be identified as a proton (highest PID "
        "probability hypothesis) with |Vz| < 20 cm")
  .note(:antineutron_emc_selection,
        "antineutron identified from the most energetic EMC shower with E > 0.48 GeV, "
        "more than 20 hit crystals, second moment S > 18 cm^2, and an opening angle "
        "greater than 20 deg with respect to any charged track")
  .note(:tag_candidate_selection,
        "when several tag candidates survive, the one with minimum |deltaE| is kept")
  .note(:efficiency_curve,
        "antineutron detection efficiency corrected data-driven as a function of the "
        "nbar momentum and cos(theta) using the J/psi -> p nbar pi- control sample "
        "at 3.097 GeV (708_3097)")
  .note(:qqbar_sideband,
        "q qbar background estimated from the MBC sideband (2.200, 2.260) GeV/c^2")
  .note(:background_veto,
        "Lambda_c+ Lambda_c- background modelled by inclusive Lambda_c- -> pbar X MC")

# Render the tag specification (takes no Selection argument) and run
alg.apply

root_files = alg.execute_on(real_data + incMC + exMCs_signal + exMCs_bkg)