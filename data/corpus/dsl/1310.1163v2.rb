# =============================================================================
# BESIII : Observation of Zc(3885) in e+e- -> pi± (D D*)∓ at sqrt(s) = 4.26 GeV
# arXiv:1310.1163v2
#
# Partial-reconstruction technique : only the bachelor pion and ONE final-state
# D meson are detected; the D* is inferred from energy-momentum conservation.
# Two parallel analyses (isospin channels) :
#   (I)  pi+ D0 D*-  with D0 -> K- pi+      (pi+ D0 tag)
#   (II) pi- D+ D*0  with D+ -> K- pi+ pi+  (pi- D+ tag)
# Data : 525 pb^-1 at 4.260 GeV
# BOSS part only : dataset preparation + event selection up to the fit.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets
### ---------------------------------------------------------------------------
data_4260  = DatasetManager.real_data.find("703_4260")     # 4.260 GeV
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")  # generic inclusive MC

# =============================================================================
# Channel I : e+e- -> pi+ D0 D*- , D0 -> K- pi+
# =============================================================================
decay_card_piD0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ D0 anti-D*-   PHSP;
  Enddecay

  Decay D0
  1.0000 K- pi+   PHSP;
  Enddecay

  End
DECAYCARD

exMC_piD0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_4260_piD0Dstar"
  config.related_dataset = data_4260
  config.events          = 100_000
  config.decay_card      = decay_card_piD0
  config.cross_section   = :default
end

alg_piD0 = Algorithm.new("Zc3885PiD0Tag")
alg_piD0.set_header(["Zc3885PiD0TagAlg/Zc3885PiD0Tag.h"])
        .set_constant({"ECMS" => [:double, 4.260]})

sel_piD0 = Selection.new
sel_piD0.select_track {
           cos_theta 0.93    # |cos(theta)| < 0.93
           Vz        10.0    # PCA < 10 cm in the beam direction
           Vr        1.0     # PCA < 1 cm in the plane perpendicular to the beam
           nTot      ">=3"   # at least three well reconstructed charged tracks
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify :kaon, against: [:pion, :proton]   # K- of the D0
           identify :pion, against: [:kaon, :proton]   # pi+ of the D0 and the bachelor pi+
           nkm  ">=1"   # at least one negatively charged track identified as a kaon
           npip ">=2"   # at least two positively charged tracks identified as pi+
         }
         # Two-constraint kinematic fit: the K- pi+ invariant mass is constrained
         # to m(D0) and the mass recoiling from the pi+ D0 system to m(D*-).
         # Only the combination with the smallest chi2 is retained.
         .kinematic_fit([:km, :pip, :pip]) {
            nominal
            invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)
            chi2_cut 200
         }

alg_piD0.note(:kinematic_fit_chi2, "Paper requires chi2 < 30 for the two-constraint fit (loose default 200 used in BOSS; the published cut is applied in the ROOT stage). The paper's fit constrains BOTH the D0 mass and the mass recoiling from pi+ D0 to m(D*-); the recoil-mass constraint on the undetected D*- is not expressible as a DSL fit constraint.")
        .note(:d0_mass_window, "K- pi+ combinations with |M(K- pi+) - m(D0)| < 15 MeV/c^2 are D0 candidates; when several combinations exist the one with mass closest to m(D0) is retained.")
        .note(:bachelor_pion_selection, "When more than one bachelor pion candidate exists, the one with the smallest chi2 of the kinematic fit is retained.")
        .note(:post_fit_cut, "M(pi+ D0) > 2.02 GeV is required to reject e+e- -> D*+ D*- , D*+ -> pi+ D0 events (a cut on the kinematically fitted four-momentum; applied in the ROOT stage).")
        .note(:recoil_mass_variable, "The measured observable is M_recoil = RM(pi D) + M(D) - m(D), the D mass resolution being minimised this way; RM is inferred from four-momentum conservation.")
        .note(:d_dbar1_background, "The D D1(2420) contribution (which would produce a near-threshold reflection) is estimated from the asymmetry of |cos(theta_piD)| and by fits to the |cos(theta_piD)|-binned data.")
        .note(:cross_feed_channel, "Cross feed from pi+ Zc(3885)- , Zc(3885)- -> D- D*0 , D*0 -> gamma/pi0 D0 (the tagged D0 being a decay product of the D*0) is simulated and included in the fit.")

alg_piD0.with_decay_card(decay_card_piD0).apply(sel_piD0)
alg_piD0.execute_on([data_4260, incMC_4260, exMC_piD0])

# =============================================================================
# Channel II : e+e- -> pi- D+ D*0 , D+ -> K- pi+ pi+
# =============================================================================
decay_card_piDp = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi- D+ anti-D*0   PHSP;
  Enddecay

  Decay D+
  1.0000 K- pi+ pi+   PHSP;
  Enddecay

  End
DECAYCARD

exMC_piDp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_4260_piDplusDstar0"
  config.related_dataset = data_4260
  config.events          = 100_000
  config.decay_card      = decay_card_piDp
  config.cross_section   = :default
end

alg_piDp = Algorithm.new("Zc3885PiDplusTag")
alg_piDp.set_header(["Zc3885PiDplusTagAlg/Zc3885PiDplusTag.h"])
        .set_constant({"ECMS" => [:double, 4.260]})

sel_piDp = Selection.new
sel_piDp.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nTot      ">=3"
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify :kaon, against: [:pion, :proton]   # K- of the D+
           identify :pion, against: [:kaon, :proton]   # two pi+ of the D+ and the bachelor pi-
           nkm  ">=1"
           npip ">=2"
           npim ">=1"   # the additional bachelor pi-
         }
         # Two-constraint kinematic fit: the K- pi+ pi+ invariant mass is constrained
         # to m(D+) and the mass recoiling from pi- D+ to m(D*0).
         .kinematic_fit([:km, :pip, :pip, :pim]) {
            nominal
            invariant_mass_of(:km, :pip, :pip).constrain_to_nominal_mass_of(:D_plus)
            chi2_cut 200
         }

alg_piDp.note(:kinematic_fit_chi2, "Paper requires chi2 < 30 for the two-constraint fit (loose default 200 used in BOSS; the published cut is applied in the ROOT stage). The paper's fit constrains BOTH the D+ mass and the mass recoiling from pi- D+ to m(D*0); the recoil-mass constraint on the undetected D*0 is not expressible as a DSL fit constraint.")
        .note(:dplus_mass_window, "K- pi+ pi+ combinations with |M(K- pi+ pi+) - m(D+)| < 15 MeV/c^2 are D+ candidates.")
        .note(:bachelor_pion_selection, "The selection is the same as for the pi+ D0 tag, with the additional requirement of a pi- track identified as the bachelor pion; when more than one bachelor pion candidate exists the one with the smallest chi2 of the kinematic fit is retained.")
        .note(:post_fit_cut, "M(pi D) requirement on the kinematically fitted four-momentum is applied in the ROOT stage.")
        .note(:cross_feed_channel, "Cross feed from pi- Zc(3885)+ , Zc(3885)+ -> D0 D*+ , D*+ -> pi0 D+ (the tagged D+ being a decay product of the D*+) is simulated and included in the fit.")
        .note(:isospin_consistency, "The two isospin channels are analysed in parallel and their sigma x B results combined by weighted average assuming isospin symmetry.")

alg_piDp.with_decay_card(decay_card_piDp).apply(sel_piDp)
alg_piDp.execute_on([data_4260, incMC_4260, exMC_piDp])
