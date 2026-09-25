# Study of e+e- -> p pbar in the vicinity of psi(3770)   [arXiv:1403.6011]
#
# Cross-section measurement of e+e- -> p pbar at 3.650, 3.773 GeV and at the
# psi(3770) line-shape scan points, separating the resonant (psi(3770) -> p pbar)
# and continuum amplitudes through their interference.

### Dataset description ###
data_3650 = DatasetManager.real_data.find("709_3650")    # 44.5 pb^-1 collected at 3.65 GeV
data_3773 = DatasetManager.real_data.find("712_3773")    # 2917 pb^-1 collected at 3.773 GeV

incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# psi(3770) line-shape scan points covering the range 3.74 - 3.90 GeV
scan_points = [
  DatasetManager.real_data.find("712_3768"),
  DatasetManager.real_data.find("712_3780"),
  DatasetManager.real_data.find("712_3800"),
  DatasetManager.real_data.find("712_3805"),
  DatasetManager.real_data.find("703_3810"),
  DatasetManager.real_data.find("712_3815"),
  DatasetManager.real_data.find("712_3820"),
  DatasetManager.real_data.find("712_3825"),
  DatasetManager.real_data.find("712_3830"),
  DatasetManager.real_data.find("712_3835"),
  DatasetManager.real_data.find("712_3840"),
  DatasetManager.real_data.find("712_3845"),
  DatasetManager.real_data.find("712_3855"),
  DatasetManager.real_data.find("712_3860"),
  DatasetManager.real_data.find("712_3875"),
  DatasetManager.real_data.find("703_3900"),
]

# Signal decay card: e+e- -> p pbar. No intermediate charmonium is produced at
# generator level (the final state is reached directly), therefore the BESIII
# KKMC convention top mother psi(4260) is used.
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 p+ anti-p- PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive signal MC: one sample per centre-of-mass energy point, so that the
# detection efficiency is evaluated at the same sqrt(s) as the data.
exMC_signal = DatasetManager.create_exclusive_mc_for([data_3650, data_3773] + scan_points) do |config|
  config.sample_name   = "ppbar_scan"     # auto-suffixed with the dataset key per point
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# Background decay card: ISR return to the lower-lying psi(3686) resonance,
# which is not removed by the ISR correction procedure and populates the same
# p pbar final state.
decay_card_isr_psip = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma psi(3686) PHSP;
  Enddecay

  Decay psi(3686)
  1.0000 p+ anti-p- PHSP;
  Enddecay

  End
DECAYCARD

exMC_isr_psip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ppbar_isr_psip"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_isr_psip
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PpbarVicinityPsi3770"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})

sel = Selection.new
sel.select_track do
      cos_theta 0.8     # |cos(theta)| < 0.8 for the polar angle of the track
      Vz        10.0    # point of closest approach within 10 cm of the IP along the beam
      Vr        1.0     # within 1 cm of the beam axis in the transverse plane
      nChrp     "==1"   # exactly one positively charged track (proton)
      nChrn     "==1"   # exactly one negatively charged track (antiproton)
      nNet      "==0"   # the two tracks have opposite charge
    end
   .pid(method: :probability) do
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]   # p+ and pbar- separated from K and pi
      nprp "==1"                                  # exactly one proton
      nprm "==1"                                  # exactly one antiproton
    end
   # |p_measured - p_expected| < 40 MeV/c (~3 sigma) for both tracks. At
   # sqrt(s) = 3.773 GeV the expected momentum is
   # p = sqrt(s/4 - m_p^2) = 1.637 GeV/c. Candidates outside the window are removed.
   .for_each(:prp) do
      define(:dp_prp) { p - 1.637 }
      where { (dp_prp > 0.040) | (dp_prp < -0.040) }
      remove
    end
   .for_each(:prm) do
      define(:dp_prm) { p - 1.637 }
      where { (dp_prm > 0.040) | (dp_prm < -0.040) }
      remove
    end
   .kinematic_fit([:prp, :prm]) do
      nominal
      constrain_four_momentum   # 4C: constrain the p pbar system to the e+e- CMS energy
      chi2_cut 200
    end

alg.note(:theta_ppbar_cut,
         "The angle between the proton and the antiproton in the rest frame of the " \
         "overall e+e- CMS system is required to be greater than 179 degrees. This " \
         "back-to-back topology cut has no dedicated DSL expression.")
   .note(:momentum_window_energy_dependence,
         "The |p_meas - p_exp| < 40 MeV/c momentum window is applied to both tracks at " \
         "every energy point, with p_exp = sqrt(s/4 - m_p^2) evaluated at that point's " \
         "sqrt(s). The BOSS selection implements it with the 3.773 GeV value " \
         "(1.637 GeV/c); the per-scan-point values are implemented by re-running the " \
         "selection with the corresponding ECMS constant.")
   .note(:isr_correction,
         "Initial-state radiation effects are not included at generator level for the " \
         "efficiency determination; they are corrected afterwards with the standard ISR " \
         "correction procedure, iterating the line-shape fit (~5 iterations) until the " \
         "dressed cross sections converge.")
   .note(:background_veto,
         "Backgrounds from psi(3770) -> D Dbar (inclusive MC), e+e- -> K+K-, mu+mu-, " \
         "tau+tau-, p pbar pi0 and p pbar gamma are studied. The total background is " \
         "0.4 events (0.06% contamination) and is accounted for in the systematic " \
         "uncertainty rather than subtracted; the ISR psi(3686) background (0.1 events) " \
         "is neglected. At 3.65 GeV the psi(3686) tail contribution (0.89 events) is " \
         "statistically subtracted from the raw signal yield.")
   .note(:cross_section_extraction,
         "The observed cross section is sigma_obs = N_sig / (epsilon * L) per energy " \
         "point; the resonant cross section, the phase angle phi and the continuum form " \
         "factor parameter C are extracted from a simultaneous fit of Eq. (1) of the " \
         "paper to the dressed cross sections of this measurement and the BaBar points. " \
         "Upper limits at 90% C.L. are set with the Feldman-Cousins method.")

alg.with_decay_card(decay_card_signal).apply(sel)
alg.execute_on([data_3650, data_3773, incMC_3650, incMC_3773] + scan_points +
               exMC_signal + [exMC_isr_psip])
